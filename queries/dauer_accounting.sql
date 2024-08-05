SELECT 
    (SELECT SUM(Dauer_culled_ind) 
     FROM dynamic_stage_transition) +
    (SELECT SUM(Dauer_to_larva) 
     FROM stage_transition) AS 'Dauer Out',

    (SELECT SUM(Dauer_culled_ind) 
     FROM dynamic_stage_transition) AS 'Dauer Culled',

    (SELECT SUM(Dauer_to_larva) 
     FROM stage_transition) AS 'Dauer to Larva',

    (SELECT Dauer_end_of_simulation_ind 
     FROM dynamic_stage_transition
     ORDER BY timestep DESC 
     LIMIT 1) AS 'Dauer Remaining';
