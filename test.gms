$title gerrymandering for whole wisconsin wards in a file

* limit printout
option mip=gurobi;
option limrow=0;
option limcol=0;
option solprint=off;


set
    wards district unit
    xy coordinate of districts
    party republicans and democrats
    center district center /d1*d8/
    ;
      
scalar
    PopUp population upper bound/1.1/
    PopLo population lower bound/0.9/
    district_num number of district(8)
    ;

district_num = card(center);
  
alias(wards,i,j,id);
alias(center,dc);

*####################### data loading starts ####################### 
*###################################################################
*###################################################################
*load voting data

$gdxin vote_data.gdx
$load wards=dim1 party=dim2

parameter vote(wards,party);

$load vote
$gdxin


*load coordinate data

$gdxin coordinate_data.gdx
$load xy=dim2
parameter coord(wards,xy);
$load coord
$gdxin

*load adjacency data
set border(i,j) yes if adjacent otherwise no;

$gdxin border_new.gdx
$load border
$gdxin

*###################################################################
*###################################################################
*####################### data loading ends #########################

*----------------------------------------------------------------

*####################### parameter setting starts ################## 
*###################################################################
*################################################################### 

* make centroid choice for each district ****to-do!!!

set centroid(i,dc) the centroid choice;
    
centroid(i,dc) = no;
centroid('F1','d1') = yes;
centroid('F1000','d2') = yes;
centroid('F2000','d3') = yes;
centroid('F3000','d4') = yes;
centroid('F4000','d5') = yes;
centroid('F5000','d6') = yes;
centroid('F6000','d7') = yes;
centroid('F7000','d8') = yes;

* create parameters storing population
parameter
    n(i) population of each ward
;
scalar    
    total_n total_population
;

n(i) = sum(party,vote(i,party));
total_n = sum(i,n(i));

* create parameters storing distance

parameter dist(i,dc) distance between ward i and center j;
dist(i,dc) = sqrt(sum(xy, sqr(coord(i,xy) - sum(j$centroid(j,dc), coord(j,xy)))));

*###################################################################
*###################################################################
*####################### parameter setting endss ###################

*----------------------------------------------------------------

*####################### variable setting starts ################### 
*###################################################################
*###################################################################

******************* variables for Model 1 *******************
free variable objectiveVal;

binary variable assign(i,dc) 1 if ward i is assigned to district center dc ;
*fix centroid point
assign.fx(i,dc)$centroid(i,dc) = 1;

nonnegative variable NV(dc) population assigned to each center;
*define population bounds
NV.UP(dc) = PopUp*sum(i,n(i))/district_num;
NV.LO(dc) = PopLo*sum(i,n(i))/district_num;

******************* extra sets and variables for Model 2 **************
set
    chosen(i) chosen border points
    to_fix(i) fixed points
    limit_dc(i,dc) limit the isolated and border points so that they can only choose the color that are adjacent to them
;


free variable
    reassignOBJ
;

positive variable
    max_pop store the maximum population of all the district
    min_pop store the minimum population of all the district
;


*###################################################################
*###################################################################
*####################### variable setting endss ####################

*----------------------------------------------------------------

*####################### model setting starts ###################### 
*###################################################################
*###################################################################

******************* Equations for Model 1 *******************
equation
    OBJ minimize population weighted distance objective function
    assgin_constr(i) each ward should be assigned to exactly one centroid ward
    population_calculation(dc) it is used to calculate the population
    ;
    
OBJ..
    objectiveVal =e= sum((i,dc), n(i)*dist(i,dc)*assign(i,dc));

assgin_constr(i)..
    sum(dc,assign(i,dc)) =e= 1;
    
population_calculation(dc)..
    NV(dc) =e= sum(i, assign(i,dc)*n(i));
    

model initialAssign /OBJ , assgin_constr , population_calculation/;

******************* Equations for Model 2 *******************

equation
    OBJ2 minimize the difference between max population and min population
    max_pop_constr(dc) calculate max population
    min_pop_constr(dc) calculate min population
    reassign(i,dc) reassignment of isolated and border points
;

OBJ2..
    reassignOBJ =e= max_pop - min_pop;

max_pop_constr(dc)..
    max_pop =g= NV(dc);

min_pop_constr(dc)..
    min_pop =L= NV(dc);

reassign(chosen,dc)..
    assign(chosen,dc) =L= sum(j$(border(chosen,j)),assign(j,dc));
    
model populationReassign /OBJ2 , assgin_constr , population_calculation, reassign, max_pop_constr, min_pop_constr/;

solve initialAssign using mip minimizing objectiveVal;

execute_unload 'initial_assignment.gdx', assign,centroid;

parameter
    isolated_reassign(i,dc) store the result of isolated points assignment by python script
    neighbor_statistic(i,dc) use this parameter to find the border point
    tempPop(dc) population for each district
;


*start dealing with isolating point and population bound

set loopNum /l1*l10/;
scalar which /0/;

loop( loopNum,
*relax the fixed points after the first loop
    assign.up(i,dc)$(not centroid(i,dc)) = 1;
    assign.lo(i,dc)$(not centroid(i,dc)) = 0;
*run script and load new assignment without isolated wards
    execute '="/Users/walter/opt/miniconda3/bin/python3" isolated.py';
    execute_load 'isolated_reassign.gdx', isolated_reassign=isolated_reassign;

    tempPop(dc) = sum(i,n(i) * isolated_reassign(i,dc));
    
    which = 0;
    break$((smax(dc,tempPop(dc)) <= PopUp*total_n/district_num) and (smin(dc,tempPop(dc)) >= PopLo*total_n/district_num));
    which = 1;
*check if a ward is on border
    neighbor_statistic(i,dc) = 0;
    neighbor_statistic(i,dc)$(sum(border(i,j),isolated_reassign(j,dc))>0.5) = 1;

*chosen = yes means the ward is on border
    chosen(i) = no;
    chosen(i)$(sum(dc,neighbor_statistic(i,dc)) > 1.5) = yes;

*fix unchosen point
    to_fix(i) = yes;
    to_fix(i) = to_fix(i) - chosen(i);
    assign.fx(to_fix(i),dc) = isolated_reassign(i,dc);

*limit the choice of border point
    limit_dc(chosen,dc) = no;
    limit_dc(chosen,dc)$(sum(j$(border(chosen,j)),isolated_reassign(j,dc))<0.5) = yes;
    assign.fx(limit_dc) = 0;

    solve populationReassign using mip minimizing reassignOBJ;
    execute_unload 'initial_assignment.gdx', assign,centroid;

);

if(which = 0,
    execute_unload 'final_assignment.gdx', isolated_reassign,centroid;
else
    execute_unload 'final_assignment.gdx', assign,centroid;
);

*display '%gams.sysdir%';
