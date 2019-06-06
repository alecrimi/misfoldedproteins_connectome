function compute_pet_prob(files)
% files is the name of the folders with files

if  nargin<1
    files = dir('.');
end

dirFlags = [files.isdir];
subFolders = files(dirFlags);

for ll = 3 : length(subFolders)
	subFolders(ll).name
	cd(subFolders(ll).name);

	sub = dir(pwd);
	dirsubFlags = [sub.isdir];
	subsub = sub(dirsubFlags);
	cd(subsub(3).name)
	sub = dir(pwd);
	dirsubFlags = [sub.isdir];
	subsub = sub(dirsubFlags);
	cd(subsub(3).name)
	sub = dir(pwd);
	dirsubFlags = [sub.isdir];
	subsub = sub(dirsubFlags);
	cd(subsub(3).name)

	pet_data = load_untouch_nii('out.nii.gz');
	pet_img = pet_data.img;

	atlas_data = load_untouch_nii('atlas_reg2.nii.gz');
	atlas_img = atlas_data.img;

	% ROIs for the Shen Atlas representing the cerebellum
	cerebellum = [118, 156,177 97 106 160 122 139 112 130 122 169 97  47 74 9 27 57  53 35  48 27 46  79 12 46 4];
	mean_val = [];
	[r,c,d] = size(atlas_img);
	for rr = 1 : r
	    for cc = 1 : c
		for dd = 1 : d

		    if( any(cerebellum == atlas_img(rr,cc,dd) )   ) 
		        mean_val(end+1) = pet_img(rr,cc,dd);
		   end

		end     
	    end
	end
	ref_val = mean(mean_val); 

	n_areas = 184; % Shen atlas has 184 ROIS
	tot_areas = zeros(n_areas,1);

	for ii = 1 : n_areas
	temp = [];
	for rr = 1 : r
	    for cc = 1 : c
		for dd = 1 : d

		    if( atlas_img(rr,cc,dd)==ii )
		       temp(end+1) = pet_img(rr,cc,dd);
		   end

		end     
	    end
	tot_areas(ii) = mean(temp);

	end

	end

	prob = tot_areas - ref_val;  % This should be improved
	prob(prob<=0)=0;
	prob
	csvwrite(['pet_' subFolders(ll).name],prob);

	cd ..
	cd ..
	cd ..
	cd ..
end
