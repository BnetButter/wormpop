SELECT
    L1Arrest_to_Larva_Mass AS 'L1 Arrest to Larva Mass',
    L1Arrest_Died_of_Starvation_Mass AS 'L1Arrest Died of Starvation Mass',
    L1Arrest_Died_of_Cull_Mass AS 'L1Arrest Died of Cull Mass',
    L1Arrest_End_of_Simulation_Mass AS 'L1 Arrest End Of Simulation Mass',
    Total_L1_Out_Mass AS 'Total L1 Out Mass',
    Total_L1_In_Simulation_Mass AS 'Total L1 In Simulation Mass',
    (L1Arrest_to_Larva_Mass * 100.0 / Total_L1_Out_Mass) AS 'L1 Arrest to Larva %',
    (L1Arrest_Died_of_Starvation_Mass * 100.0 / Total_L1_Out_Mass) AS 'L1Arrest Died of Starvation %',
    (L1Arrest_Died_of_Cull_Mass * 100.0 / Total_L1_Out_Mass) AS 'L1Arrest Died of Cull %',
    (Total_L1_Out_Mass * 100.0 / Total_L1_In_Simulation_Mass) AS 'Total L1 Out Mass %',
    (L1Arrest_End_of_Simulation_Mass * 100.0 / Total_L1_In_Simulation_Mass) AS 'End Of Simulation Mass %'
FROM (
    SELECT
        (SELECT SUM(l1arrest_to_larva_mass) FROM stage_transition) AS L1Arrest_to_Larva_Mass,
        (SELECT SUM(L1Arrest_starvation_mass) FROM dynamic_stage_transition) AS L1Arrest_Died_of_Starvation_Mass,
        (SELECT SUM(L1Arrest_culled_mass) FROM dynamic_stage_transition) AS L1Arrest_Died_of_Cull_Mass,
        (SELECT SUM(L1Arrest_end_of_simulation_mass) FROM dynamic_stage_transition) AS L1Arrest_End_of_Simulation_Mass,
        (
            (SELECT SUM(l1arrest_to_larva_mass) FROM stage_transition) +
            (SELECT SUM(L1Arrest_starvation_mass) FROM dynamic_stage_transition) +
            (SELECT SUM(L1Arrest_culled_mass) FROM dynamic_stage_transition)
        ) AS Total_L1_Out_Mass,
        (
            SELECT SUM(Max_Mass)
            FROM (
                SELECT Worm_Name, MAX(Mass) AS Max_Mass 
                FROM worms 
                WHERE stage LIKE '%L1%'
                GROUP BY Worm_Name
            ) AS MaxMassTable
        ) AS Total_L1_In_Simulation_Mass
) AS Subquery;

