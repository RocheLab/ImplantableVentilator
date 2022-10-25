function [Output] = analyzeSnipPFV(FlowData,VolAutoData,Vollocs,PActlocs,Flowlocs,Vminbounds)
%Analyse Snip for peak flows, peak volumes, and the distance between peaks
%of flows, volumes and actuator pressure


% initialize the first row of the output matrix
Output = zeros(1,10); % content is [bound1 bound2 Fpkloc Fpkval Ppkloc Vpkloc Vpkval F_P P_V F_V];


% use vol min as the first set of bounds to check

bound1 = Vminbounds(1);%initialize
bound2 = Vminbounds(2);%initialize
row = 1; %initialize

boundlookup = Vminbounds; %initialize

%start with a set of volume minimum bounds

while bound1 < boundlookup(end) %may need to change this as the for loop bounds
    Vol_i = find(Vollocs>bound1&Vollocs<bound2);%find the Vol peak in between two indices of volume bounds
    PAct_i = find(PActlocs>bound1&PActlocs<bound2);%find the PAct peak in between two indices of volume bounds
    Flow_i = find(Flowlocs>bound1&Flowlocs<bound2);% find the Flow pk
    if isempty(Vol_i) && isempty(PAct_i) && isempty(Flow_i)
        % skip this loop
        % reassign bounds
        nextbound = find(boundlookup>bound2);
        bound1 = bound2;
        bound2 = boundlookup(nextbound(1));
        continue
    elseif length(Vol_i) == 1 && length(PAct_i_i) == 1 && length(Flow_i) == 1
        % assume this set out bounds is fine, assign appropriate values
        Fpkloc = Flowlocs(Flow_i);
        Fpkval = FlowData(Fpkloc);
        Ppkloc = PActlocs(PAct_i);
        Vpkloc = Vollocs(Vol_i);
        Vpkval = VolAutoData(Vpkloc);
        
        F_P = Fpkloc - Ppkloc; % pos if PAct peak is first, neg if Flow peak is first
        P_V = Ppkloc - Vpkloc; % pos if Vol peak is first, neg if PAct peak is first
        F_V = Fpkloc - Vpkloc; % pos if Vol peak is first, neg if Flow peak is first
        
        Output(row,:) = [bound1 bound2 Fpkloc Fpkval Ppkloc Vpkloc Vpkval F_P P_V F_V];
        %reassign for next for loop
        row = row+1;
        nextbound = find(boundlookup>bound2);
        bound1 = bound2;
        bound2 = boundlookup(nextbound(1));
        continue
    elseif isempty(Vol_i) == 1 %if there is no peak volume between the two indices, assume even if there is vol, it's too low to detect
        % skip this loop
        % reassign bounds
        nextbound = find(boundlookup>bound2);
        bound1 = bound2;
        bound2 = boundlookup(nextbound(1));
        continue
    elseif  isempty(PAct_i) == 1 &&  isempty(Flow_i) == 1 %probably false max
        % skip this loop
        % reassign bounds
        nextbound = find(boundlookup>bound2);
        bound1 = bound2;
        bound2 = boundlookup(nextbound(1));
        continue
    elseif length(Vol_i) == 1 && (length(Flow_i)>1 || length(PAct_i)>1)% either PAct of F act has more than one but only 1 vol pk
        %take the closest pact or flow pk
        
        Vpkloc = Vollocs(Vol_i);
        Vpkval = VolAutoData(Vpkloc);
        if isempty(Flow_i)
             [Fpkval,FMidx] = max(FlowData(bound1:Vpkloc)); %find a max flow before the Vpk loc
             Fpkloc = FMidx+bound1-1; 
        elseif length(Flow_i)==1
            Fpkloc = Flowlocs(Flow_i);
            Fpkval = FlowData(Fpkloc);
        else %many flow pks
              %find closest flow pk before vol max
              listF = find(Flowlocs(Flow_i)<Vpkloc);
              Fpkloc = Flowlocs(Flow_i(listF(end)));
              Fpkval = FlowData(Fpkloc);
        end
        
        if isempty(PAct_i)
            Ppkloc = 0; %no pact Pk
            F_P = NaN; % pos if PAct peak is first, neg if Flow peak is first
            P_V = NaN; % pos if Vol peak is first, neg if PAct peak is first
        elseif length(PAct_i) == 1
            Ppkloc = PActlocs(PAct_i);
            
            F_P = Fpkloc - Ppkloc; % pos if PAct peak is first, neg if Flow peak is first
            P_V = Ppkloc - Vpkloc; % pos if Vol peak is first, neg if PAct peak is first
        else %many Pacts
            %take Pact closes to Vpk
            [[],Pind] = min(abs(PActlocs(PAct_i)-Vpkloc));
            Ppkloc = PActlocs(PAct_i(Pind));

            F_P = Fpkloc - Ppkloc; % pos if PAct peak is first, neg if Flow peak is first
            P_V = Ppkloc - Vpkloc; % pos if Vol peak is first, neg if PAct peak is first
        end

        F_V = Fpkloc - Vpkloc; % pos if Vol peak is first, neg if Flow peak is first
        
        Output(row,:) = [bound1 bound2 Fpkloc Fpkval Ppkloc Vpkloc Vpkval F_P P_V F_V];
        %reassign for next for loop
        row = row+1;
        nextbound = find(boundlookup>bound2);
        bound1 = bound2;
        bound2 = boundlookup(nextbound(1));
        continue
        
        
    elseif length(Vol_i)>1 
        % narrow bounds
        newbound =0.5*Vollocs(Vol_i(2))-0.5*Vollocs(Vol_i(1)); %set new bound 2 athalfwaypoint between vpks
        if newbound > bound2
        %skip
        else
            b1ind = find(boundlookup==bound1);
            boundlookup = [boundlookup(1:b1ind) newbound boundlookup(b1ind+1:end)]; %insert new bound into boundlookup
            %bound1 = bound1;
            bound2 = newbound;
        end
        continue
  
   
        
    elseif length(Vol_i)==1 && (length(Flow_i)==1 || length(PAct_i)==1) % if there is 1 vol peak and either 1 F or 1 P pk and 0 of the other
        if length(Flow_i)==1 && isempty(PAct_i)% a spontaneous breath

            Fpkloc = Flowlocs(Flow_i);
            Fpkval = FlowData(Fpkloc);
            Ppkloc = 0; %no pact Pk
            Vpkloc = Vollocs(Vol_i);
            Vpkval = VolAutoData(Vpkloc);
            F_P = NaN; %no pact Pk
            P_V = NaN ;%no pact Pk
            F_V = Fpkloc - Vpkloc; % pos if Vol peak is first, neg if Flow peak is first

            Output(row,:) = [bound1 bound2 Fpkloc Fpkval Ppkloc Vpkloc Vpkval F_P P_V F_V];
            %reassign for next for loop
            row = row+1;
            nextbound = find(boundlookup>bound2);
            bound1 = bound2;
            bound2 = boundlookup(nextbound(1));
            continue
        elseif length(PAct_i) == 1 && isempty(Flow_i)% a very small flow peak that went below detection (prob destructive breath) 
                % assume this set of bounds is fine, assign appropriate values
                Vpkloc = Vollocs(Vol_i);
                [Fpkval,FMidx] = max(FlowData(bound1:Vpkloc)); %find a max flow before the Vpk loc
                Fpkloc = FMidx+bound1-1; 
                Ppkloc = PActlocs(PAct_i);

                Vpkval = VolAutoData(Vpkloc);

                F_P = Fpkloc - Ppkloc; % pos if PAct peak is first, neg if Flow peak is first
                P_V = Ppkloc - Vpkloc; % pos if Vol peak is first, neg if PAct peak is first
                F_V = Fpkloc - Vpkloc; % pos if Vol peak is first, neg if Flow peak is first

                Output(row,:) = [bound1 bound2 Fpkloc Fpkval Ppkloc Vpkloc Vpkval F_P P_V F_V];
                %reassign for next for loop
                row = row+1;
                nextbound = find(boundlookup>bound2);
                bound1 = bound2;
                bound2 = boundlookup(nextbound(1));
                continue
        
        
        end
    end
end
end

