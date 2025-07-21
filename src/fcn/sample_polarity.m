function [polarity, vote] = sample_polarity(Igray)
% Assesses the sample polarity in an image
% Imput: a gray-scale image
% The sample image is divided into sub image
% Output: String: polarity, 'bright' or 'dark'. 
%         Integer: vote, 0 means unanimously dark, 9 means unaniously
%         bright
% Author: Oscar Sten, Istituto Italiano di Technologia, 13/12/2021.
[n, m] = size(Igray);
T_list = [];
mean_color_list = [];
for i=1:3
    for j=1:3
        sub_img = Igray(1+floor(n*(i-1)/3):floor((n*i)/3), ...
            (1+floor(m*(j-1)/3):floor((m*j)/3)));
        [T, ~] = otsuthresh(imhist(sub_img));
        mean_color = mean(sub_img(:))/255;
        T_list = [T_list; T];
        mean_color_list = [mean_color_list; mean_color];
    end
end
polarity_votes = (mean_color_list < T_list);
vote = sum(polarity_votes);
if vote < 4.5
    polarity = 'dark';
else
    polarity = 'bright';
end

end