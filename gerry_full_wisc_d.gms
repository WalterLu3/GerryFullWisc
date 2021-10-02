$title gerrymandering for whole wisconsin wards

option mip=gurobi;
option limrow=0;
option limcol=0;
option solprint=off;


set wards
    xy
    party
    center district center /d1*d8/
    ;
      
scalar
    PopUp population upper bound/1.1/
    PopLo population lower bound/0.9/
    district_num
    ;
    
district_num = card(center);
  
alias(wards,i,j,id);
alias(center,dc);
*load voting data
$gdxin vote_data.gdx
$load wards=dim1 party=dim2

parameter vote(wards,party);
$load vote
$gdxin

parameter n(i) population of each ward;

n(i) = sum(party,vote(i,party));

*load coordinate data

$gdxin coordinate_data.gdx
$load xy=dim2
parameter coord(wards,xy);
$load coord
$gdxin


set
    centroid(i,dc) the centroid choice
    contig(i,dc) yes if contiguous constraint is needed 
    ;
    
    
centroid(i,dc) = no;
centroid('F1','d1') = yes;
centroid('F1000','d2') = yes;
centroid('F2000','d3') = yes;
centroid('F3000','d4') = yes;
centroid('F4000','d5') = yes;
centroid('F5000','d6') = yes;
centroid('F6000','d7') = yes;
centroid('F7000','d8') = yes;

contig(i,dc) = yes;
contig(centroid) = no;



parameter dist(i,dc) distance between ward i and center j;
dist(i,dc) = sqrt(sum(xy, sqr(coord(i,xy) - sum(j$centroid(j,dc), coord(j,xy)))));


*parameter distance(i,j) distance between ward i and j;

*load distance data
*$gdxin distance.gdx
*$load distance=d
*$gdxin
*$offtext
$ontext
set	closer(i,j,dc)	Flag that ward i is closer to distict dc centroid than is ward j;
$gdxin closer.gdx
$load closer
$gdxin
$offtext


*load adjacency data
set border(i,j) yes if adjacent otherwise no;

*parameter adj(i,j);

$gdxin border.gdx
$load border
$gdxin

loop(centroid(i,dc), contig(j,dc)$border(j,i) = no;);
nonnegative variable
    NV(dc) population assigned to each center;

binary variable
    assign(i,dc) 1 if ward i is assigned to district center dc ;
    
assign.fx(i,dc)$centroid(i,dc) = 1;

free variable
    objectiveVal;
    

    
   
equation
    OBJ population weighted distance objective function
    assgin_constr(i) each ward should be assigned to exactly one centroid ward
    population_calculation(dc)
*    contiguous_constr(i,dc)
    ;
    
OBJ..
    objectiveVal =e= sum((i,dc), n(i)*dist(i,dc)*assign(i,dc));

assgin_constr(i)..
    sum(dc,assign(i,dc)) =e= 1;
    
population_calculation(dc)..
    NV(dc) =e= sum(i, assign(i,dc)*n(i));
    
*contiguous_constr(i,dc)$(contig(i,dc))..
*    assign(i,dc) =L= sum(j$(border(i,j) and closer(j,i,dc)),assign(j,dc));
    

*contiguity



model wardHess /all/;

NV.UP(dc) = PopUp*sum(i,n(i))/district_num;
NV.LO(dc) = PopLo*sum(i,n(i))/district_num;

solve wardHess using mip minimizing objectiveVal;


set
    to_assign(i,dc);
    
to_assign(i,dc)$(assign.l(i,dc) >= 0.5) = yes;

display to_assign;

execute_unload 'fullWisc_Wo_Adj_constraint_weighted_objective.gdx', assign;
    
 