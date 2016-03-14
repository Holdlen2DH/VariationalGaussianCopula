function vec=post2vec(post)
% convert posterior stucture to a 1=-d col. vector
vec=[post.m;post.c(:)];
