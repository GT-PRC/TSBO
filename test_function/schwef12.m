%Schwefel1.2
% Function Description:
%      xmin  = [0, 0, 0.....0]  (all zeroes)
%      fxmin = 0                  (zero)
function Schw = schwef12(Swarm)
[SwarmSize, Dim] = size(Swarm);
% Schw=zeros(SwarmSize,1);
Schw=(Swarm(:,1:1)).^2;
for i=2:Dim
Schw = Schw+((sum(Swarm(:,1:i)')).^2)';
end
%-100 100
%-100 100