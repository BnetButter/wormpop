-- COUNT number of L1

SELECT
    L1Arrest_to_Larva AS 'L1 Arrest to Larva',
    L1Arrest_Died_of_Starvation AS 'L1Arrest Died of Starvation',
    L1Arrest_Died_of_Cull AS 'L1Arrest Died of Cull',
    L1Arrest_End_of_Simulation AS 'L1 Arrest End Of Simulation',
    Total_L1_Out AS 'Total L1 Out',
    (L1Arrest_to_Larva * 100.0 / Total_L1_Out) AS 'L1 Arrest to Larva %',
    (L1Arrest_Died_of_Starvation * 100.0 / Total_L1_Out) AS 'L1Arrest Died of Starvation %',
    (L1Arrest_Died_of_Cull * 100.0 / Total_L1_Out) AS 'L1Arrest Died of Cull %',
	(L1Arrest_End_of_Simulation / Total_L1_In_Simulation * 100) AS 'L1 At End of Simulation %'
FROM (
    SELECT
	    (SELECT COUNT(DISTINCT Worm_Name) FROM worms WHERE stage LIKE '%L1%') AS Total_L1_In_Simulation,
        (SELECT SUM(L1Arrest_to_Larva) FROM stage_transition) AS L1Arrest_to_Larva,
        (SELECT SUM(L1Arrest_starvation_ind) FROM dynamic_stage_transition) AS L1Arrest_Died_of_Starvation,
        (SELECT SUM(L1Arrest_culled_ind) FROM dynamic_stage_transition) AS L1Arrest_Died_of_Cull,
        (SELECT SUM(L1Arrest_end_of_simulation_ind) FROM dynamic_stage_transition) AS L1Arrest_End_of_Simulation,
        (
            (SELECT SUM(L1Arrest_to_Larva) FROM stage_transition) +
            (SELECT SUM(L1Arrest_starvation_ind) FROM dynamic_stage_transition) +
            (SELECT SUM(L1Arrest_culled_ind) FROM dynamic_stage_transition)
        ) AS Total_L1_Out
) AS Subquery;
