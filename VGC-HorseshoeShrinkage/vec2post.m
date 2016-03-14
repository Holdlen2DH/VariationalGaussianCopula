function post=vec2post(vec,postDim)
% convert vector of posterior params to posterior structure
post.m=vec(1:postDim);
vec=vec(postDim+1:end);
post.c=reshape(vec(:),postDim,postDim);
end
