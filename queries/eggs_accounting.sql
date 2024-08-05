SELECT 
    (SELECT SUM(Total_Eggs_Laid) 
     FROM worm_summary) AS Total_Eggs_Laid,

    (SELECT SUM(Total_Eggs_Laid) * 65 
     FROM worm_summary) AS Total_Eggs_Laid_MASS,

    (SELECT COUNT(*) 
     FROM worms 
     WHERE stage = 'dead' 
     AND Notes LIKE '%cull%' 
     AND Mass <= 65) AS Eggs_Removed_By_Culling,

    (SELECT COUNT(*) * 65 
     FROM worms 
     WHERE stage = 'dead' 
     AND Notes LIKE '%cull%' 
     AND Mass <= 65) AS Eggs_Removed_By_Culling_MASS,

    (SELECT COUNT(*) 
     FROM worm_summary 
     WHERE Larva_span_days IS NOT NULL) AS Hatched,

    (SELECT COUNT(*) * 65 
     FROM worm_summary 
     WHERE Larva_span_days IS NOT NULL) AS Hatched_MASS
