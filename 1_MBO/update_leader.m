function [leader,Z_leader,followings,Z_following]=update_leader(leader,Z_leader,followings,Z_following)
followings=[followings;leader];
Z_following=[Z_following,Z_leader];
leader=followings(1,:);
Z_leader=Z_following(1);
followings(1,:)=[];
Z_following(1)=[];