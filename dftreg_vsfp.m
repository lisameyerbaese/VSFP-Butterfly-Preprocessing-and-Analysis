function [A_reg, D_reg, reg_info] = dftreg_vsfp(A_input, D_input, upfac, ref_imageA, ref_imageD)

[sX,sY,sZ] = size(A_input);
A_reg = zeros(sX,sY,sZ);
D_reg = zeros(sX,sY,sZ);

[reg_info, fft_reg] = dftregistration(fft2(ref_imageA), fft2(ref_imageD), upfac);
ref_imageD_reg = abs(ifft2(fft_reg));

for i = 1:sZ
    [~,fft_reg] = dftregistration(fft2(ref_imageA), fft2(A_input(:,:,i)), upfac);
    A_reg(:,:,i) = abs(ifft2(fft_reg));
end

for i = 1:sZ
    [~,fft_reg] = dftregistration(fft2(ref_imageD), fft2(D_input(:,:,i)), upfac);
    D_reg(:,:,i) = abs(ifft2(fft_reg));
end
end
