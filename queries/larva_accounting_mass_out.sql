SELECT
    Larva_to_Adult AS 'Larva to Adult ng',
    Larva_Culled AS 'Larva Culled ng',
    Larva_to_Dauer AS 'Larva to Dauer ng',
    Larva_starved AS 'Larva starved ng',
    Larva_to_inefficiency AS 'Larva to inefficiency ng',
    Larva_to_Metabolic_Cost AS 'Larva to Metabolic Cost ng',
    Larva_remaining AS 'Larva remaining ng',
    Mass_Larva_Out AS 'Mass Larva Out ng',
    (Larva_to_Adult * 100.0 / Mass_Larva_Out) AS 'Larva to Adult %',
    (Larva_Culled * 100.0 / Mass_Larva_Out) AS 'Larva Culled %',
    (Larva_to_Dauer * 100.0 / Mass_Larva_Out) AS 'Larva to Dauer %',
    (Larva_starved * 100.0 / Mass_Larva_Out) AS 'Larva starved %',
    (Larva_to_inefficiency * 100.0 / Mass_Larva_Out) AS 'Larva to inefficiency %',
    (Larva_to_Metabolic_Cost * 100.0 / Mass_Larva_Out) AS 'Larva to Metabolic Cost %',
    (Larva_remaining * 100.0 / Mass_Larva_Out) AS 'Larva remaining %'
FROM (
    SELECT
        (SELECT SUM(larva_to_adult_mass) FROM stage_transition) AS Larva_to_Adult,
        (SELECT SUM(Larva_culled_mass) FROM dynamic_stage_transition) AS Larva_Culled,
        (SELECT SUM(larva_to_dauer_mass) FROM stage_transition) AS Larva_to_Dauer,
        (SELECT SUM(larva_starvation_mass) FROM dynamic_stage_transition) AS Larva_starved,
        (SELECT SUM(Metabolic_Efficiency_Loss) FROM worms WHERE stage = 'larva') AS Larva_to_inefficiency,
        (SELECT SUM(Metabolic_Cost) FROM worms WHERE stage = 'larva') AS Larva_to_Metabolic_Cost,
        (SELECT SUM(Larva_end_of_simulation_mass) FROM dynamic_stage_transition) AS Larva_remaining,
        (
            (SELECT SUM(larva_to_adult_mass) FROM stage_transition) +
            (SELECT SUM(Larva_culled_mass) FROM dynamic_stage_transition) +
            (SELECT SUM(larva_to_dauer_mass) FROM stage_transition) +
            (SELECT SUM(larva_starvation_mass) FROM dynamic_stage_transition) +
            (SELECT SUM(Metabolic_Efficiency_Loss) FROM worms WHERE stage = 'larva') +
            (SELECT SUM(Metabolic_Cost) FROM worms WHERE stage = 'larva') +
            (SELECT SUM(Larva_end_of_simulation_mass) FROM dynamic_stage_transition)
        ) AS Mass_Larva_Out
) AS Subquery
