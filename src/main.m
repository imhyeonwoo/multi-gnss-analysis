%% GNSS Experiment - Minimal Pipeline (Step 1)
% 목적: 한 세션 데이터를 LS로 위치추정 → ECEF→LLH/ENU → 기본 플롯/KML
% 디렉토리 구조(루트 기준):
%  - data/Team3_Novatel_0918.mat, Team3_Android_0918_01.mat
%  - src/ (이 파일과 LS.m, xyz2llh.m, xyz2enu.m, eleazi.m, skyplot_sys.m)
%  - out/ (자동 생성)

clc; clear; close all;

%% [A] 경로/세션 설정
here      = fileparts(mfilename('fullpath'));  % src
root_dir  = fileparts(here);                   % 루트
data_dir  = fullfile(root_dir, 'data');
out_dir   = fullfile(root_dir, 'out');
if ~exist(out_dir,'dir'); mkdir(out_dir); end
addpath(here);  % src를 path에 추가

% --- 사용할 세션 선택(필요시 아래 하나만 주석 해제) ---
% dataset = 'android';
dataset = 'novatel';

switch lower(dataset)
    case 'novatel'
        data_file = fullfile(data_dir, 'Team3_Novatel_0918.mat');
        kml_name  = 'novatel_0918_baseline.kml';
    case 'android'
        data_file = fullfile(data_dir, 'Team3_Android_0918_01.mat');
        kml_name  = 'android_0918_01_baseline.kml';
    otherwise
        error('dataset must be "novatel" or "android"');
end

assert(exist('LS.m','file')==2, 'LS.m이 src에 있어야 합니다.');
assert(isfile(data_file), '데이터 파일을 찾을 수 없습니다: %s', data_file);

% 기준점 LLH [rad, rad, m]; 모르면 NaN(평균 위치를 기준으로 사용)
cfg.ref_llh = [NaN NaN NaN];

%% [B] 데이터 로드(유연 로더: 변수명 자동 감지)
obs = loadObsMat(data_file);
% obs 구조: t [Nx1], pr [NxM], sat [1xM], sys [1xM], (옵션) el/az/cn0
fprintf('[B] Loaded %s | epochs=%d | sats=%d\n', ...
        dataset, numel(obs.t), numel(obs.sat));

%% [C] LS 위치추정
N  = numel(obs.t);
Ms = numel(obs.sat);
x_ecef = nan(N,3); clk = nan(N,1);
min_sats = 4;

for k = 1:N
    idx = true(1,Ms);                  % 1단계: 전체 위성 사용(마스크 없음)
    if sum(idx) < min_sats || any(isnan(obs.pr(k,idx))); continue; end
    range_k = obs.pr(k,idx);
    sat_k   = obs.sat(idx);
    sys_k   = obs.sys(idx);
    [x_hat, clk_bias] = LS(range_k, sat_k, sys_k);  % 수업 제공 LS.m 시그니처 가정
    x_ecef(k,:) = x_hat(:).';
    clk(k)      = clk_bias;
end

valid = all(isfinite(x_ecef),2);
assert(any(valid), '유효 추정 epoch이 없습니다. 데이터/LS 입력을 확인하세요.');
fprintf('[C] LS done. valid=%d/%d\n', sum(valid), N);

%% [D] 좌표 변환(ECEF→LLH/ENU)
llh = nan(N,3); enu = nan(N,3);

% 기준점 결정(미지정 시 유효 추정 평균 사용)
if any(isnan(cfg.ref_llh))
    ref_llh = xyz2llh(mean(x_ecef(valid,:),1).');
else
    ref_llh = cfg.ref_llh(:);
end
ref_xyz = llh2xyz_local(ref_llh);

for k = find(valid).'
    llh(k,:) = xyz2llh(x_ecef(k,:).').';
    enu(k,:) = xyz2enu(x_ecef(k,:).', ref_xyz).';
end
fprintf('[D] Transform done.\n');

%% [E] ENU 플롯 저장
fig = figure('Name', sprintf('ENU_%s', dataset), 'Color','w');
t   = obs.t(valid);
subplot(3,1,1); plot(t, enu(valid,1), '.-'); ylabel('E (m)'); grid on;
subplot(3,1,2); plot(t, enu(valid,2), '.-'); ylabel('N (m)'); grid on;
subplot(3,1,3); plot(t, enu(valid,3), '.-'); ylabel('U (m)'); xlabel('epoch'); grid on;
saveas(fig, fullfile(out_dir, sprintf('%s_enu.png', dataset)));

%% [F] KML 저장(구글어스 확인용)
writeKML_simple(fullfile(out_dir, kml_name), llh(valid,:));
fprintf('[F] KML written: %s\n', fullfile(out_dir, kml_name));

%% [G] 요약
fprintf(['[Summary]\n',...
    '- Dataset: %s\n- Valid epochs: %d\n- Outputs: ENU plot, KML\n'], ...
    dataset, sum(valid));

%% ======== 로컬 함수들 ========

function obs = loadObsMat(mat_path)
% 다양한 .mat 형식을 견디는 유연 로더(필수 누락 시 안전한 기본값 생성)
    S = load(mat_path);

    % 0) 최상위에 구조체 하나만 있는 경우 그걸 펼치기
    if numel(fieldnames(S))==1
        fn = fieldnames(S); maybe = S.(fn{1});
        if isstruct(maybe) || istable(maybe)
            S = structify(maybe);
        end
    end

    F = struct();
    % 필수 1: pseudorange 행렬
    F.pr  = pickVar(S, {'pr','pseudorange','rho','range','RHO','PR'}, false);
    pr = F.pr;
    if isempty(pr)
        error('필수 변수(pr/rho/range)가 필요합니다.');
    end
    % 방향 정리: [N x M]
    if size(pr,1) < size(pr,2), pr = pr.'; end
    [N,M] = size(pr);

    % 필수 2: 시간 t (없으면 1:N로 자동 생성)
    t = pickVar(S, {'t','time','epoch','gpsTime','GPST','T'}, true);
    if isempty(t), t = (1:N).'; end
    t = t(:);
    if numel(t) ~= N
        % 차원 안 맞으면 앞 N개만 사용
        t = t(1:N);
    end

    % 권장: sat/sys (없으면 임시 생성)
    sat = pickVar(S, {'sat','svid','prn','sv','sat_id','PRN'}, true);
    if isempty(sat), sat = 1:M; end
    sat = sat(:).';

    sys = pickVar(S, {'sys','system','gnss','const','constellation'}, true);
    if isempty(sys), sys = ones(1,M); end
    sys = sys(:).';

    % 옵션: el/az/cn0
    el  = pickVar(S, {'el','elev','elevation'}, true);
    az  = pickVar(S, {'az','azim','azimuth'}, true);
    cn0 = pickVar(S, {'cn0','CNo','cnr','snr','SNR'}, true);

    obs = struct('t',t, 'pr',pr, 'sat',sat, 'sys',sys);
    if ~isempty(el),  obs.el  = shapeOpt(el,N,M);  end
    if ~isempty(az),  obs.az  = shapeOpt(az,N,M);  end
    if ~isempty(cn0), obs.cn0 = shapeOpt(cn0,N,M); end

    % ---- 내부 보조 ----
    function X = shapeOpt(X,N,M)
        if isvector(X)
            X = repmat(X(:).', N, 1);  % [1xM]이면 N번 복제
        end
        if size(X,1)~=N || size(X,2)~=M
            % 가능한 한 맞춰보기
            if size(X,1)==M && size(X,2)==N, X = X.'; end
            X = X(1:min(end,N), 1:min(end,M));
            if size(X,1)<N || size(X,2)<M
                X(N,M) = NaN;  % 부족분 NaN 채움
            end
        end
    end

    function S2 = structify(x)
        % table/timetable → struct, struct는 그대로
        if istable(x)
            S2 = struct();
            for vn = string(x.Properties.VariableNames)
                S2.(vn) = x.(vn);
            end
        else
            S2 = x;
        end
    end
end

function v = pickVar(S, names, optional)
% names 중 첫 번째로 존재하는 변수를 반환. 없으면 [](optional) 또는 에러.
    if nargin<3, optional=false; end
    v = [];
    % 최상위 필드 검사
    for i=1:numel(names)
        if isfield(S, names{i}), v = S.(names{i}); return; end
    end
    % 구조체 안쪽 1단계 탐색
    fns = fieldnames(S);
    for f = 1:numel(fns)
        val = S.(fns{f});
        if isstruct(val)
            for i=1:numel(names)
                if isfield(val, names{i}), v = val.(names{i}); return; end
            end
        elseif istable(val)
            for i=1:numel(names)
                if any(strcmp(val.Properties.VariableNames, names{i}))
                    v = val.(names{i});
                    return;
                end
            end
        end
    end
    if ~optional
        error('Variable not found: one of %s', strjoin(names, ', '));
    end
end

