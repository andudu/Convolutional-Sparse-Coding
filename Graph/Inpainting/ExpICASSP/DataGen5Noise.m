% generate and save relavant data for inpainting test

% load the image
S0 = imresize(double(stdimage('lena.grey'))/255,.5);

noise_level = [.4,.5,.6,.65,.7,.75];   %5 different noise levels

for i = 1: length(noise_level)
    
    perc_noise = noise_level(i) ;
    %generate missing pixels
    ind = [];
    ind(1,:) = randi(size(S0,1),1,ceil(size(S0,1)*size(S0,2)*perc_noise));
    ind(2,:) = randi(size(S0,2),1,ceil(size(S0,1)*size(S0,2)*perc_noise));
    S_c = S0; 
    for t = 1:size(ind,2)
        S_c(ind(1,t),ind(2,t)) = 0;
    end
    save(['CorImNoise256',num2str(i),'.mat'],'S_c','ind');
    
end

