% Necessite un add on 
function [] = dispFeature(A)
sz = size(A);
act1 = reshape(A,[sz(1) sz(2) 1 sz(3)]);
I = imtile(mat2gray(act1), 'GridSize',[8 12]);
imshow(I)
end