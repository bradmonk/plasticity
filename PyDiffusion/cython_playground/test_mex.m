zap = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
disp('zap before:');
disp(zap);
zap_rows = 3;
zap_new = foo_mex(zap, zap_rows);
disp('foo_mex(zap, 3):');
disp(zap_new);
disp('zap after:');
disp(zap);
