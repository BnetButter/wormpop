SELECT 
    (SELECT SUM(Dauer_culled_ind) FROM dynamic_stage_transition) +
    (SELECT SUM(Dauer_to_larva) FROM stage_transition) AS 'Dauer Out';



SELECT SUM(Dauer_culled_ind) as 'Dauer Culled' FROM dynamic_stage_transition;
SELECT SUM(Dauer_to_larva) as 'Dauer to Larva' FROM stage_transition;

SELECT Dauer_end_of_simulation_ind as 'Dauer Remaining' FROM dynamic_stage_transition
ORDER BY timestep DESC LIMIT 1;

