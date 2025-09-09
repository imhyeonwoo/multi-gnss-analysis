function [ ECEF, CLK ] = LS( range, sat, sys )
len_SVN = length(range); % 위치계산 사용 위성 수
sys_unique = unique(sys); % 위치계산 사용 위성의 Satellite System 추출
if len_SVN >= length(sys_unique)+3 % 시스템의 수에 따라 최소 필요 위성수 만족하는 경우 위치 계산 수행
    num_sys = length(sys_unique); % 시스템 수
    % 초기화
    state = zeros(3+num_sys,1); % 처음에 수신기의 위치를 알수 없으므로 임의의 값으로 지정
    X = ones(3+num_sys,1); % State 초기화
    
    % Satelite System이 4개 모두 들어오지 않을 경우 H 행렬과 선형화 지점 정리하기 위한 sys_mat 행렬 선연
    if num_sys ~= 4
        sys_mat = nan(len_SVN,1);
        for n = 1:num_sys
            sys_mat(sys == sys_unique(n)) = n;
        end
    else
        sys_mat = sys;
    end
    
    count = 1;
    % 오차범위가 0.001보다 작아질때까지 최소자승법을 반복 사용
    while(norm(X) > 0.001)
        count = count + 1;
         if count > 500
             ECEF = nan(1,3);
              CLK = nan(1,4);
            break
         end
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 필수 내용 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H = ;       % H 관측 행렬 생성
        Y = ;       % Y 행렬 생성
        X = ;       % 최소자승법 수행
        state = ;   % State 업데이트
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    ECEF = state(1:3)'; % ECEF 좌표 계산 결과 저장
    CLK = nan(1,4); % 시계 오차 계산 결과 변수 선언
    CLK(sys_unique) = state(4:end)'; % 시계 오차 계산 결과 저장
else % % 시스템의 수에 따라 최소 필요 위성수를 만족하지 못한 경우 nan 저장
    ECEF = nan(1,3);
    CLK = nan(1,4);
end


