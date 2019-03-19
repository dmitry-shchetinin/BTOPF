function [A, b] = BTenvelope_for_sqrt_xy( indices, varbounds, Vdif, dims )
%returns sparse matrix and structure with vectors of lower and upper bounds
%for linear envelopes of z=sqrt(xy) function. Each envelope is represented as
%z=a1*x+a2*y+b, l<b<u

%unpack data
lb=(varbounds.lb).^2;
ub=(varbounds.ub).^2;
ind_x=indices.x;
ind_y=indices.y;
L=dims.nrows;
%initialize rhs vector and vectors for constructing sparse matrix
ind_vals=zeros(2*L,1);
b=struct('l', zeros(L,1), 'u', zeros(L,1));
for i=1:L
    %obtain coefficients
    [coefs, offsets]=compute_sqrt_xy_bounds(lb(ind_x(i)),ub(ind_x(i)),lb(ind_y(i)),ub(ind_y(i)), Vdif.min(i), Vdif.max(i));
    %record results
    %x
    ind_vals(i)=coefs(1);
    %y
    ind_vals(L+i)=coefs(2);
    %offsets
    b.l(i)=offsets.min;
    b.u(i)=offsets.max;
end
A=sparse(repmat(transpose(1:L),2,1), [ind_x; ind_y], ind_vals, L, dims.ncols);
end


%compute quadratic upper bound and linear lower bound for sqrt(xy)
function [coefs, offsets]=compute_sqrt_xy_bounds(xl, xu, yl, yu, offset_min, offset_max)
%get intersections of two boundary lines
point1=intersections_of_line_with_box(offset_min, xl, xu, yl, yu);
point2=intersections_of_line_with_box(offset_max, xl, xu, yl, yu);
%get intersections of middle line between boundary lines or box corners with box
if (point1.x1~=point2.x1 && point1.y1~=point2.y1)
    point3.x1=xl;
    point3.y1=yl;
else
    point3.x1=(point1.x1+point2.x1)/2;
    point3.y1=(point1.y1+point2.y1)/2;
end
if (point1.x2~=point2.x2 && point1.y2~=point2.y2)
    point3.x2=xu;
    point3.y2=yu;
else
    point3.x2=(point1.x2+point2.x2)/2;
    point3.y2=(point1.y2+point2.y2)/2;
end

%get unique intersection points of three lines with Vbox
x_temp1=[]; y_temp1=[]; x_temp2=[]; y_temp2=[];
if (point1.x1~=point1.x2 || point1.y1~=point1.y2)
    x_temp1=point1.x2; y_temp1=point1.y2;
end
if (point2.x1~=point2.x2 || point2.y1~=point2.y2)
    x_temp2=point2.x2; y_temp2=point2.y2;
end
xall=[point1.x1; point2.x1; point3.x1; point3.x2; x_temp1; x_temp2];
yall=[point1.y1; point2.y1; point3.y1; point3.y2; y_temp1; y_temp2];
zall=sqrt(xall.*yall);

%compute parameters of bounds that have form z=a1*x+a2*y+b
x_temp=(point3.x1+point3.x2)/2;
y_temp=(point3.y1+point3.y2)/2;
coefs=[sqrt(y_temp/x_temp)/2; sqrt(x_temp/y_temp)/2];
offsets.max=sqrt(x_temp*y_temp)-coefs(1)*x_temp-coefs(2)*y_temp;
offsets.min=min(zall-coefs(1)*xall-coefs(2)*yall);
end


function point=intersections_of_line_with_box(offset, xl, xu, yl, yu)
    %beginning of the line segment
    if (offset<=xl-yl)
        point.x1 = xl;
        point.y1 = xl - offset;
    else
        point.x1 = yl + offset;
        point.y1 = yl;
    end
    %end of the line segment
    if (offset>=xu-yu)
        point.x2 = xu;
        point.y2 = xu - offset;
    else
        point.x2 = yu + offset;
        point.y2 = yu;
    end
end