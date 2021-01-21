function [asc, ord]=areapiena(a,b,dim);
    asc1=[1:dim];
    asc2=[dim:-1:1];
    asc=[asc1 asc2];
    for i=1:dim
        ord(i)=a(i);
        ord(dim+i)=b((dim+1)-i);
    end

end