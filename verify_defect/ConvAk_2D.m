%Computes the k-ic discrete convolution A*A*...*A
function convAk = ConvAk_2D(A, k)
	%Extract the truncation order
	M = length(A)-1;
	
	%Pad each dimension of A to the next power of 2 of the minimum padding (k-1)M+1 in the (+,+) quadrant
	len_pad = 2^(nextpow2(2*k*M+1)-1) - (M+1);
	A_pad = padarray(A, [len_pad, len_pad], 0, 'post');
	
	%Mirror the (+,+) quadrant to the full Fourier coefficients
	len_pad = length(A_pad);
	A_pad = padarray(padarray(A_pad, [len_pad, len_pad], 'symmetric', 'pre'), [1, 1], 0, 'pre');
	A_pad(:, len_pad+1)=[]; A_pad(len_pad+1, :)=[];
	
	if ~isa(A, 'intval')
		%Without interval arithmetic, we can directly use Matlab's fft and ifft functions
		convAk = ifft2(fft2(ifftshift(A_pad)).^k);
	else
		%With interval arithmetic, verifyfft must be used instead
		convAk = verifyfft(verifyfft(ifftshift(A_pad), 1).^k, -1);
	end
	
	%Only return the first (1+kM)^2 non-zero components
	convAk = 1.0 * convAk(1:1+k*M, 1:1+k*M);
end
