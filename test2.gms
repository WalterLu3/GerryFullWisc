$title gerrymandering for whole wisconsin wards in a file

* limit printout
option mip=gurobi;
option limrow=0;
option limcol=0;
option solprint=off;


scalar
    starttime
    total_run_time;

starttime = jnow;


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
parameter order(i) used to randomly chose the district centers;
order(i) = 1;
embeddedCode python:
    import random
    order = list(gams.get("order"))
    values = list(range(1,7079))
    random.shuffle(values)
    new_order = []
    for idx in range(0,7078):
        new_order.append((order[idx][0],values[idx]))
    gams.set("order",new_order)
endEmbeddedCode order


set centroid(i,dc) the centroid choice
    original_centroid(i,dc) ;
  
  
centroid(i,dc) = no;
*centroid(i,dc)$((order(i) < (ord(dc)+0.5)) and (order(i) > (ord(dc)-0.5))) = yes;


$ontext
centroid('F5665','d1') = yes;
centroid('F5390','d2') = yes;
centroid('F2015','d3') = yes;
centroid('F6603','d4') = yes;
centroid('F5359','d5') = yes;
centroid('F2140','d6') = yes;
centroid('F2296','d7') = yes;
centroid('F845','d8') = yes;
$offtext


centroid('F6026','d1') = yes;
centroid('F5298','d2') = yes;
centroid('F4764','d3') = yes;
centroid('F4285','d4') = yes;
centroid('F5097','d5') = yes;
centroid('F1726','d6') = yes;
centroid('F6076','d7') = yes;
centroid('F6687','d8') = yes;

$ontext
centroid("F2463","d1") = yes;
centroid("F5926","d2") = yes;
centroid("F5902","d3") = yes;
centroid("F1777","d4") = yes;
centroid("F5734","d5") = yes;
centroid("F2411","d6") = yes;
centroid("F2887","d7") = yes;
centroid("F7005","d8") = yes;
$offtext
original_centroid(i,dc) = centroid(i,dc); 
* create parameters storing population
parameter
    n(i) population of each ward
;
scalar    
    total_n total_population
    LowerBound pop lower bound
    UpperBound pop Upper bound
;

n(i) = sum(party,vote(i,party));
total_n = sum(i,n(i));

UpperBound = PopUp*total_n/district_num;
LowerBound = PopLo*total_n/district_num;


* create parameters storing distance

parameter dist_i_j(i,j) distance between ward i and j
          dist(i,dc) distance between ward i and center dc;;

$gdxin distance_between_i_j.gdx
$load dist_i_j=dist
$gdxin


dist(i,dc) = sum(j$centroid(j,dc),dist_i_j(i,j));
*parameter dist(i,dc) distance between ward i and center j;
*dist(i,dc) = sqrt(sum(xy, sqr(coord(i,xy) - sum(j$centroid(j,dc), coord(j,xy)))));


*###################################################################
*###################################################################
*####################### parameter setting endss ###################

*----------------------------------------------------------------

*##################center change model starts#######################
*###################################################################
*###################################################################

set ward_subset(i) the subset of wards
    temp_centroid(i,dc);

ward_subset(i) = no;

binary variable middle(i) the binary variable to choose new district center;

free variable objDist;

equations
    unique_middle
    objCenterChange
;

unique_middle..
    sum(i$ward_subset(i),middle(i)) =e= 1;
objCenterChange..
    sum((i,j)$(ward_subset(i) and ward_subset(j)),middle(i)*dist_i_j(i,j)) =e= objDist;
    

model changeCentroid /unique_middle, objCenterChange/;

*solve changeCentroid using mip minimizing objDist;

*set middle_point(i);

*middle_point(i)$(middle.l(i) > 0.5) = yes;

*display middle_point;


*###################################################################
*###################################################################
*##################center change model ends#########################

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
NV.UP(dc) = UpperBound;
NV.LO(dc) = LowerBound;

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
    centroid_neighbor(i,dc) at least one districgt is adjacent to centroid
    ;
    
OBJ..
    objectiveVal =e= sum((i,dc),n(i)*dist(i,dc)*assign(i,dc));

assgin_constr(i)..
    sum(dc,assign(i,dc)) =e= 1;
    
population_calculation(dc)..
    NV(dc) =e= sum(i, assign(i,dc)*n(i));
    
*centroid_neighbor(i,dc)$centroid(i,dc)..
*    sum(j$border(i,j), assign(j,dc)) =g= sum(j$border(i,j), 1);

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

*relax population bound
NV.UP(dc) = +inf;
NV.LO(dc) = 0;

parameter
    isolated_reassign(i,dc) store the result of isolated points assignment by python script
    neighbor_statistic(i,dc) use this parameter to find the border point
    tempPop(dc) population for each district
;

scalar upup, lolo, upas, loas;

*start dealing with isolating point and population bound

set loopNum /l1*l5/;
parameter inner_border(i,dc) see if a ward is on innder border of a district
          outer_border(i,dc) see if a ward is on outer border of a district
          average_border_distance(dc) average border distance form dc to its chosen border points
          largest_border_distance(dc)
;
scalar which /0/;

File popFile /exp_change_centroid.csv/;

popFile.ap = 1;

loop( loopNum,
    temp_centroid(i,dc) = no;
*change the centroid to a new ward
    loop( dc,
        ward_subset(i) = no;
*get the wards of a district
        ward_subset(i)$(assign.l(i,dc) > 0.5) = yes;
        solve changeCentroid using mip minimizing objDist;
        temp_centroid(i,dc)$(middle.l(i)>0.5 and assign.l(i,dc) > 0.5) = yes;
    );
    centroid(i,dc) = temp_centroid(i,dc);
    execute_unload 'initial_assignment.gdx', assign,centroid;

*relax the fixed points after the first loop
    assign.up(i,dc) = 1;
    assign.lo(i,dc) = 0;
*run script and load new assignment without isolated wards
    execute '="/Users/walter/opt/miniconda3/bin/python3" isolated.py';
    execute_load 'isolated_reassign.gdx', isolated_reassign=isolated_reassign;

    tempPop(dc) = sum(i,n(i) * isolated_reassign(i,dc));
$ontext    
    upup = PopUp*total_n/district_num;
    lolo = PopLo*total_n/district_num;
    upas = smax(dc,tempPop(dc));
    loas = smin(dc,tempPop(dc));
    display upup,lolo,upas,loas;
$offtext
    which = 0;
    break$((smax(dc,tempPop(dc)) <= PopUp*total_n/district_num) and (smin(dc,tempPop(dc)) >= PopLo*total_n/district_num));

    which = 1;
*check if a ward is on border
    neighbor_statistic(i,dc) = 0;
    neighbor_statistic(i,dc)$(sum(border(i,j),isolated_reassign(j,dc))>0.5) = 1;

*check if a ward is on inner or outer border
    inner_border(i,dc) = 0;
    outer_border(i,dc) = 0;
*    inner_border(i,dc)$(sum(center,neighbor_statistic(i,center)) > 1.5 and isolated_reassign(i,dc)>0.5
*                        and tempPop(dc) > PopUp*total_n/district_num) = 1;
*    outer_border(i,dc)$(neighbor_statistic(i,dc) > 0.5 and isolated_reassign(i,dc)<0.5
*                        and tempPop(dc) < PopLo*total_n/district_num) = 1;

    inner_border(i,dc)$(sum(center,neighbor_statistic(i,center)) > 1.5 and isolated_reassign(i,dc)>0.5
                        and tempPop(dc) > PopUp*total_n/district_num ) = 1;
    outer_border(i,dc)$(neighbor_statistic(i,dc) > 0.5 and isolated_reassign(i,dc)<0.5
                        and tempPop(dc) < PopLo*total_n/district_num ) = 1;
    
*    inner_border(i,dc)$(sum(center,neighbor_statistic(i,center)) > 1.5 and isolated_reassign(i,dc)>0.5
*                        and abs(tempPop(dc)-UpperBound) < 0.1*(UpperBound -LowerBound) ) = 1;
*    outer_border(i,dc)$(neighbor_statistic(i,dc) > 0.5 and isolated_reassign(i,dc)<0.5
*                        and abs(tempPop(dc)-LowerBound) < 0.1*(UpperBound -LowerBound)) = 1;
    
*    break$(ord(loopNum) = 2);

*chosen = yes means the ward is on border
    chosen(i) = no;
*    chosen(i)$((sum(dc,neighbor_statistic(i,dc)) > 1.5) and not(sum(dc$centroid(i,dc),1) > 0.5)) = yes;
    chosen(i)$(sum(dc,inner_border(i,dc)) > 0.5) = yes;
    chosen(i)$(sum(dc,outer_border(i,dc)) > 0.5) = yes;

*calculate a dictrict border's average distance    
    average_border_distance(dc) = 0;
    average_border_distance(dc)$(sum(i$(chosen(i) and isolated_reassign(i,dc)> 0.5),1) > 0.5) = sum(i$(chosen(i) and isolated_reassign(i,dc)> 0.5),dist(i,dc))/sum(i$(chosen(i) and isolated_reassign(i,dc)> 0.5),1);
    
    largest_border_distance(dc) = 0;
    largest_border_distance(dc) = smax(i$(chosen(i) and isolated_reassign(i,dc)> 0.5),dist(i,dc));

    chosen(i)$(sum(dc$(isolated_reassign(i,dc)> 0.5 and dist(i,dc) < average_border_distance(dc)),1) > 0) = no;
*    chosen(i)$(sum(dc$(isolated_reassign(i,dc)> 0.5 and dist(i,dc) < (average_border_distance(dc)+largest_border_distance(dc))/2),1) > 0) = no;
*fix unchosen point
    to_fix(i) = yes;
    to_fix(i) = to_fix(i) - chosen(i);
    assign.fx(to_fix(i),dc) = isolated_reassign(i,dc);

*limit the choice of border point
    limit_dc(i,dc) = no;
    limit_dc(chosen,dc)$(sum(j$(border(chosen,j)),isolated_reassign(j,dc))<0.5) = yes;
    assign.fx(limit_dc) = 0;

    solve populationReassign using mip minimizing reassignOBJ;
    
*    break$(populationReassign.objVal > LowerBound);

    execute_unload 'initial_assignment.gdx', assign,centroid;
    
);

display centroid;

total_run_time = (jnow - starttime)*24*3600;

if(which = 0,
    execute_unload 'final_assignment.gdx', isolated_reassign,centroid;
    put popFile populationReassign.objVal, total_run_time, "Y";
else
    execute_unload 'final_assignment.gdx', assign,centroid;
    put popFile populationReassign.objVal, total_run_time, "N";
    execute_unload 'failed_case.gdx', original_centroid;
    display "failed";
);

*display '%gams.sysdir%';
