function imgout = poisson_blend(im_s, mask_s, im_t)
% -----Input
% im_s     source image (object)
% mask_s   mask for source image (1 meaning inside the selected region)
% im_t     target image (background)
% -----Output
% imgout   the blended image

[imh, imw, nb] = size(im_s);

% imgin = im_t .* (1-mask_s) + im_s .* mask_s;

if nb > 1
    for i = 1:nb
        imgout(:,:,i) = poisson_blend(im_s(:,:,i), mask_s, im_t(:,:,i));
    end
    return;
end

V = zeros(imh, imw);
V(1:imh*imw) = 1:imh*imw;

b = zeros(imh*imw, 1);
e = 1;
index = 1;
% total # of i, j, v = (imh-1)*(imw-1)*4 + (imh-1)*2 + (imw-1)*2 + 4; 
i = zeros(1, (imh-1)*(imw-1)*4 + (imh-1)*2 + (imw-1)*2 + 4);
j = zeros(1, (imh-1)*(imw-1)*4 + (imh-1)*2 + (imw-1)*2 + 4);
v = zeros(1, (imh-1)*(imw-1)*4 + (imh-1)*2 + (imw-1)*2 + 4);


%TODO: fill the elements in A and b, for each pixel in the image
% inner
for y = 2 : imh-1
    for x = 2 : imw-1
        i(index : index+4) = [e,      e,          e           e,          e];
        j(index : index+4) = [V(y,x), V(y,x+1),   V(y,x-1),   V(y+1,x),   V(y-1,x)];
        v(index : index+4) = [4,      -1,         -1,         -1,         -1];
        b(e, 1) = 4*im_s(y,x) - im_s(y,x+1) - im_s(y,x-1) - im_s(y+1,x) - im_s(y-1,x);

        if mask_s(y,x) == 1
            if mask_s(y,x+1) == 0
                j(index+1) = 1;
                v(index+1) = 0;
                b(e,1) = b(e,1) + im_t(y,x+1);
            end
            if mask_s(y,x-1) == 0
                j(index+2) = 1;
                v(index+2) = 0;
                b(e,1) = b(e,1) + im_t(y,x-1);
            end
            if mask_s(y+1,x) == 0
                j(index+3) = 1;
                v(index+3) = 0;
                b(e,1) = b(e,1) + im_t(y+1,x);
            end
            if mask_s(y-1,x) == 0
                j(index+4) = 1;
                v(index+4) = 0;
                b(e,1) = b(e,1) + im_t(y-1,x);
            end
        end

        index = index + 5;
        e = e + 1;
    end
end

% first row + last row:
for x = 2 : imw-1
    i(index : index+2) = [e,      e,        e];
    j(index : index+2) = [V(1,x), V(1,x+1), V(1,x-1)];
    v(index : index+2) = [2,      -1,       -1];
    b(e, 1) = 2*im_s(1,x) - im_s(1,x+1) - im_s(1,x-1);

    if mask_s(1,x) == 1
        if mask_s(1,x+1) == 0
            j(index+1) = 1;
            v(index+1) = 0;
            b(e,1) = b(e,1) + im_t(1,x+1);
        end
        if mask_s(1,x-1) == 0
            j(index+2) = 1;
            v(index+2) = 0;
            b(e,1) = b(e,1) + im_t(1,x-1);
        end
    end
    index = index + 3;
    e = e + 1;

    i(index : index+2) = [e,        e,          e];
    j(index : index+2) = [V(imh,x), V(imh,x+1), V(imh,x-1)];
    v(index : index+2) = [2,        -1,         -1];
    b(e, 1) = 2*im_s(imh,x) - im_s(imh,x+1) - im_s(imh,x-1);
    if mask_s(imh,x) == 1
        if mask_s(imh,x+1) == 0
            j(index+1) = 1;
            v(index+1) = 0;
            b(e,1) = b(e,1) + im_t(imh,x+1);
        end
        if mask_s(imh,x-1) == 0
            j(index+2) = 1;
            v(index+2) = 0;
            b(e,1) = b(e,1) + im_t(imh,x-1);
        end
    end
    index = index + 3;
    e = e + 1;
end
% first col + last col
for y = 2 : imh-1
    i(index : index+2) = [e,      e,        e];
    j(index : index+2) = [V(y,1), V(y+1,1), V(y-1,1)];
    v(index : index+2) = [2,      -1,       -1];
    b(e, 1) = 2*im_s(y,1)  - im_s(y+1,1) - im_s(y-1,1);
    if mask_s(y,1) == 1
        if mask_s(y+1,1) == 0
            j(index+1) = 1;
            v(index+1) = 0;
            b(e,1) = b(e,1) + im_t(y+1,1);
        end
        if mask_s(y-1,1) == 0
            j(index+2) = 1;
            v(index+2) = 0;
            b(e,1) = b(e,1) + im_t(y-1,1);
        end
    end
    index = index + 3;
    e = e + 1;

    i(index : index+2) = [e,        e,          e];
    j(index : index+2) = [V(y,imw), V(y+1,imw), V(y-1,imw)];
    v(index : index+2) = [2,        -1,         -1];
    b(e, 1) = 2*im_s(y,imw) - im_s(y+1,imw) - im_s(y-1,imw);

    if mask_s(y,imw) == 1
        if mask_s(y+1,imw) == 0
            j(index+1) = 1;
            v(index+1) = 0;
            b(e,1) = b(e,1) + im_t(y+1,imw);
        end
        if mask_s(y-1,imw) == 0
            j(index+2) = 1;
            v(index+2) = 0;
            b(e,1) = b(e,1) + im_t(y-1,imw);
        end
    end
    index = index + 3;
    e = e + 1;
end

i(index : index+3) = [e,      e+1,      e+2,      e+3];
j(index : index+3) = [V(1,1), V(1,imw), V(imh,1), V(imh,imw)];
v(index : index+3) = [1,      1,        1,        1];
b(e, 1) = im_s(1,1); % left top
b(e+1, 1) = im_s(1,imw); % right top
b(e+2, 1) = im_s(imh,1); % left bottom
b(e+3, 1) = im_s(imh,imw); % right bottom
e = e + 4;


assignin('base','i',i);
assignin('base','j',j);
assignin('base','v',v);

A = sparse(i,j,v);
assignin('base','A',A);
assignin('base','b',b);


%use "lscov" or "\",
solution = A\b;
error = sum(abs(A*solution-b));
disp(error)
assignin('base','solution',solution);

%      in the output image to obtain the blended image
% masked = mask_s;
% masked(masked > 0) = solution;
% imgout = im_t * (1-mask_s) + masked;
% imshow(imgout)

% f = figure(), hold off,
imgout = reshape(solution,[imh,imw]);
imgout = im_t .* (1-mask_s) + imgout .* mask_s;
% imshow(imgout)