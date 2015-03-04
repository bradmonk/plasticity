function [] = bradftp(fname)

bradftp = ftp('ftp.bradleymonk.com','user','pass');
cd(bradftp);
cd(bradftp,'bradleymonk.com/matlab')
	
	% Change path to directory with fname.m
	wfiname = which(fname);
	wfoname = regexprep(wfiname, '(\w*\.m$)', '','once');

	w1 = what(wfoname);
	pathNow = cd(w1.path);

mput(bradftp, fname);

	% Change path back to previous directory
	cd(pathNow)

close(bradftp);


%{
bradftp = ftp('ftp.bradleymonk.com','monakhos','8r4d13y!');
cd(bradftp);
cd(bradftp,'bradleymonk.com/matlab')
	
	% Change path to directory with MAINBOX.m
	w1 = what('mainfunctions');
	pathMain = w1.path;
	w2 = what('miscfunctions');
	pathMisc = w2.path;
	pathNow = cd(pathMain);

mput(bradftp, 'MAINBOX.m');

	% Change path back to previous directory
	cd(pathNow)

close(bradftp);
%}
%{
bradftp = ftp('ftp.matlabs.net','monakhos','8r4d13y!');
cd(bradftp)
cd(bradftp,'matlabs.net')
	
	% Change path to directory with MAINBOX.m
	w1 = what('mainfunctions');
	pathMain = w1.path;
	w2 = what('miscfunctions');
	pathMisc = w2.path;
	pathNow = cd(pathMain);

mput(bradftp, 'MAINBOX.m');

	% Change path back to previous directory
	cd(pathNow)

close(bradftp);

%s = regexprep(wfname, '(\.m)', '$`','once')
w1 = what('mainfunctions');
	pathMain = w1.path;
	w2 = what('miscfunctions');
	pathMisc = w2.path;
%}

end