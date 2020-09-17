function pso = PSOnew(m,d0,rho,h,r,f,d)
    %% Problem Definiton

    nVar = length(m);        % Number of Unknown (Decision) Variables

    VarSize = [1 nVar];      % Matrix Size of Decision Variables

    VarMin = 0;	             % Lower Bound of Decision Variables
    VarMax = 200;            % Upper Bound of Decision Variables


    %% Parameters of PSO

    kappa = 1;
    phi1 = 2.05;
    phi2 = 2.05;
    phi = phi1+phi2;
    chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));
    
    MaxIt = 10;        % Maximum Number of Iterations
    nPop = 1000;        % Population Size (Swarm Size)

    w = chi;             % Intertia Coefficient
    wdamp = 1;      % Damping Ratio of Inertia Coefficient
    c1 = chi*phi1;            % Personal Acceleration Coefficient
    c2 = chi*phi2;            % Social Acceleration Coefficient

    % The Flag for Showing Iteration Information
       

    MaxVelocity = 0.2*(VarMax-VarMin);
    MinVelocity = -MaxVelocity;
    
    %% Initialization

    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Cost = inf;

    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
        
        %Generate dp from random solution
        dp=forward_HEM(particle(i).Position(1:length(rho)),particle(i).Position(length(rho)+1:length(m)),h,r,f);

        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        particle(i).Cost =  sqrt(norm((d0-dp)./d0)^2/length(d0));

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end

    end

    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt, 1);


    %% Main Loop of PSO

    for it=1:MaxIt

        for i=1:nPop

            % Update Velocity
            particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
            particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
            
            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
            
            % Apply Lower and Upper Bound Limits
            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);

            %Update dp from new position
            dp=forward_HEM(particle(i).Position(1:length(rho)),particle(i).Position(length(rho)+1:length(m)),h,r,f);
            
            % Evaluation
            particle(i).Cost = sqrt(norm((d0-dp)./d0)^2/length(d0));

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end            

            end

        end

        % Store the Best Cost Value
        BestCosts(it) = GlobalBest.Cost;
        
        % Display the Best Costs
        disp(['Iteration ' num2str(it) ':Best Cost= ' num2str(BestCosts(it))]);
        
        % Damping Inertia Coefficient
        w = w * wdamp;

    end
    pso = GlobalBest.Position;
    
end

    