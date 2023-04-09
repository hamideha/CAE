clc; clear;

%% TRAINING
training_img = double(imread("training.jpg"));

r = training_img(:,:,1);
g = training_img(:,:,2);
b = training_img(:,:,3);

r_padded = padarray(r, [2 2], 'symmetric', 'both');
g_padded = padarray(g, [2 2], 'symmetric', 'both');
b_padded = padarray(b, [2 2], 'symmetric', 'both');

r = r(:);
g = g(:);
b = b(:);

Xr = im2col(r_padded, [5,5]);
Xg = im2col(g_padded, [5,5]);
Xb = im2col(b_padded, [5,5]);

[Ag_rggb, Ab_rggb] = rggb(Xr, Xg, Xb, g, b);
[Ar_gbrg, Ab_gbrg] = gbrg(Xr, Xg, Xb, r, b);
[Ar_grbg, Ab_grbg] = grbg(Xr, Xg, Xb, r, b);
[Ar_bggr, Ag_bggr] = bggr(Xr, Xg, Xb, r, g);

%% RAW MOSAIC GENERATION
% IF img_i isn't the ground truth RGB you can't calculate RMSE
img_i = imread("test_image.jpg");
img = im2double(img_i);
[m, n, z] = size(img);

% USE THIS IF DIRECTLY GIVEN BAYER PATTERN
% mosaic = im2double(imread("test1.png"));
% USE THIS IF SIMULATING A BAYER PATTERN FROM GROUND TRUTH
mosaic = bayer(img, "rggb");
imwrite(mosaic, "bayer.png");

bayer_img = padarray(mosaic, [3 3], 'symmetric', 'both');
bayer_img(end - 2, :) = []; bayer_img(:, end - 2) = [];
bayer_img(3, :) = []; bayer_img(:, 3) = [];

%% LINEAR REGRESSION
% Assume RGGB Bayer filter
final_r = zeros(m, n); 
final_r(1:2:end, 1:2:end) = mosaic(1:2:end, 1:2:end);
final_g = zeros(m, n); 
final_g(1:2:end, 2:2:end) = mosaic(1:2:end, 2:2:end);
final_g(2:2:end, 1:2:end) = mosaic(2:2:end, 1:2:end);
final_b = zeros(m, n); 
final_b(2:2:end, 2:2:end) = mosaic(2:2:end, 2:2:end);

for row = 3:m + 2
    for col = 3:n + 2
        img_mono = bayer_img(row - 2:row + 2, col - 2:col + 2);
        img_mono = img_mono(:);
        if mod(row, 2) && mod(col, 2)
            final_b(row - 2, col - 2) = Ab_rggb'*img_mono;
            final_g(row - 2, col - 2) = Ag_rggb'*img_mono;
        elseif ~mod(row, 2) && ~mod(col, 2)
            final_g(row - 2, col - 2) = Ag_bggr'*img_mono;
            final_r(row - 2, col - 2) = Ar_bggr'*img_mono;
        elseif mod(row, 2) && ~mod(col, 2)
            final_r(row - 2, col - 2) = Ar_grbg'*img_mono;
            final_b(row - 2, col - 2) = Ab_grbg'*img_mono;
        elseif ~mod(row, 2) && mod(col, 2)
            final_b(row - 2, col - 2) = Ab_gbrg'*img_mono;
            final_r(row - 2, col - 2) = Ar_gbrg'*img_mono;
        end
    end
end

final = cat(3, final_r, final_g, final_b);
final = uint8(final.*255);
native_demosaic = demosaic(imread("bayer.png"), "rggb");
imwrite(native_demosaic, "matlabdemosaic.png")
imwrite(final, "linreg.png")

figure;
imshow(final);
title("Demosaiced image using linear regression");

figure;
imshow(native_demosaic);
title("Demosaiced image using MATLAB demosaic()");

rmse = calculate_rmse(final, img_i);
fprintf("The RMSE between the original and demosaiced image is: %.5f\n", rmse);

rmse_de = calculate_rmse(native_demosaic, img_i);
fprintf("The RMSE between the original and demosaiced image using MATLAB demosaic() is: %.5f\n", rmse_de);

error = rgb2gray(final - img_i).^2;
imwrite(error, "error.png");

function rmse = calculate_rmse(img, original_img)
    rmse = sqrt(immse(img, original_img));
end

function bayer_img = bayer(img, pattern)
    [h, w, ~] = size(img);
    % Green
    bayer_img = img(:,:,2);

    if pattern == "rggb"    
        % Red
        bayer_img(1:2:h, 1:2:w) = img(1:2:h, 1:2:w, 1);
        % Blue
        bayer_img(2:2:h, 2:2:w) = img(2:2:h, 2:2:w, 3); 
    elseif pattern == "bggr"
        % Red
        bayer_img(2:2:h, 2:2:w) = img(2:2:h, 2:2:w, 1);
        % Blue
        bayer_img(1:2:h, 1:2:w) = img(1:2:h, 1:2:w, 3); 
    elseif pattern == "gbrg"
        % Red
        bayer_img(2:2:h, 1:2:w) = img(2:2:h, 1:2:w, 1);
        % Blue
        bayer_img(1:2:h, 2:2:w) = img(1:2:h, 2:2:w, 3);
    elseif pattern == "grbg"
        % Red
        bayer_img(1:2:h, 2:2:w) = img(1:2:h, 2:2:w, 1);
        % Blue
        bayer_img(2:2:h, 1:2:w) = img(2:2:h, 1:2:w, 3);
    else
        error("Invalid pattern. Must be either 'rggb', 'bggr', 'grbg', or 'gbgr'")
    end
end

function [Ag,Ab] = rggb(Xr, Xg, Xb, green, blue)
    r = zeros(5); r(1:2:end, 1:2:end) = 1;
    g = zeros(5); g(2:2:end) = 1;
    b = zeros(5); b(2:2:end, 2:2:end) = 1;

    X = (r(:).*Xr + g(:).*Xg + b(:).*Xb)';

    Ag = (X'*X)\X'*green;
    Ab = (X'*X)\X'*blue;
end

function [Ar,Ab] = gbrg(Xr, Xg, Xb, red, blue)
    r = zeros(5); r(2:2:end, 1:2:end) = 1;
    g = zeros(5); g(1:2:end) = 1;
    b = zeros(5); b(1:2:end, 2:2:end) = 1;

    X = (r(:).*Xr + g(:).*Xg + b(:).*Xb)';
    
    Ar = (X'*X)\X'*red;
    Ab = (X'*X)\X'*blue;
end

function [Ar,Ab] = grbg(Xr, Xg, Xb, red, blue)
    r = zeros(5); r(1:2:end, 2:2:end) = 1;
    g = zeros(5); g(1:2:end) = 1;
    b = zeros(5); b(2:2:end, 1:2:end) = 1;

    X = (r(:).*Xr + g(:).*Xg + b(:).*Xb)';
    
    Ar = (X'*X)\X'*red;
    Ab = (X'*X)\X'*blue;
end

function [Ar,Ag] = bggr(Xr, Xg, Xb, red, green)
    r = zeros(5); r(2:2:end, 2:2:end) = 1;
    g = zeros(5); g(2:2:end) = 1;
    b = zeros(5); b(1:2:end, 1:2:end) = 1;

    X = (r(:).*Xr + g(:).*Xg + b(:).*Xb)';
    
    Ar = (X'*X)\X'*red;
    Ag = (X'*X)\X'*green;
end
