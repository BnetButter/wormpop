SELECT
    Bacteria_into_Larva AS 'Bacteria into Larva ng',
    Eggs_to_Larva AS 'Eggs to Larva ng',
    L1_to_Larva AS 'L1 to Larva ng',
    Dauer_to_Larva AS 'Dauer to Larva ng',
    Mass_into_Larva,
    (Bacteria_into_Larva * 100.0 / Mass_into_Larva) AS 'Bacteria into Larva %',
    (Eggs_to_Larva * 100.0 / Mass_into_Larva) AS 'Eggs to Larva %',
    (L1_to_Larva * 100.0 / Mass_into_Larva) AS 'L1 to Larva %',
    (Dauer_to_Larva * 100.0 / Mass_into_Larva) AS 'Dauer to Larva %'
FROM (
    SELECT
        (SELECT SUM(Amount_Eaten) FROM worms WHERE stage = 'larva') AS Bacteria_into_Larva,
        (SELECT SUM(egg_to_larva_mass) FROM stage_transition) AS Eggs_to_Larva,
        (SELECT SUM(l1arrest_to_larva_mass) FROM stage_transition) AS L1_to_Larva,
        (SELECT SUM(dauer_to_larva_mass) FROM stage_transition) AS Dauer_to_Larva
) AS Subquery,
(
    SELECT 
        (Subquery.Bacteria_into_Larva + Subquery.Eggs_to_Larva + Subquery.L1_to_Larva + Subquery.Dauer_to_Larva) AS Mass_into_Larva
    FROM (
        SELECT
            (SELECT SUM(Amount_Eaten) FROM worms WHERE stage = 'larva') AS Bacteria_into_Larva,
            (SELECT SUM(egg_to_larva_mass) FROM stage_transition) AS Eggs_to_Larva,
            (SELECT SUM(l1arrest_to_larva_mass) FROM stage_transition) AS L1_to_Larva,
            (SELECT SUM(dauer_to_larva_mass) FROM stage_transition) AS Dauer_to_Larva
    ) AS Subquery
) AS Total
