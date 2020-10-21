function id = genID

id = [num2str(datenum(date)) '_'];
idrand = ceil(5.*rand(5,1));
for i=1:5
	id = [id num2str(idrand(i))];
end
