SELECT
    -- Count unique number of worms
    (SELECT COUNT(DISTINCT Worm_Name) FROM worm_summary) AS unique_worm_count,

    -- Count total body mass of worms in ng
    (SELECT SUM(Total_Body_Mass) FROM worm_summary) AS total_body_mass_sum,

    -- Count born as dauer
    (SELECT COUNT(*) FROM worm_summary WHERE Dauer_span_days IS NOT NULL AND Larva_span_days IS NULL) AS BORN_DAUER,

    -- Count mass of born as dauer
    (SELECT SUM(Total_Body_Mass) FROM worm_summary WHERE Dauer_span_days IS NOT NULL AND Larva_span_days IS NULL) AS BORN_DAUER_MASS_NG,

    -- Count born as egg
    (SELECT COUNT(*) FROM worm_summary WHERE Larva_span_days IS NOT NULL) AS BORN_EGG,

    -- Count mass of worms born as egg
    (SELECT SUM(Total_Body_Mass) FROM worm_summary WHERE Larva_span_days IS NOT NULL) AS BORN_EGG_MASS,

    -- Count all dead worms
    (SELECT COUNT(*) FROM worms WHERE Notes LIKE "%Cause of Death%") AS DEAD_WORM_COUNT,

    -- Count all worms died by culled
    (SELECT COUNT(*) FROM worms WHERE Notes LIKE "%culled%") AS DEAD_WORM_CULLED_COUNT,

    -- Count mass of worms died by culled
    (SELECT SUM(Mass) FROM worms WHERE Notes LIKE "%culled%") AS DEAD_WORM_CULLED_MASS,

    -- Count all worms died by starvation
    (SELECT COUNT(*) FROM worms WHERE Notes LIKE "%starve%") AS DEAD_WORM_STARVE_COUNT,

    -- Count mass of worms died by starvation
    (SELECT SUM(Mass) FROM worms WHERE Notes LIKE "%starve%") AS DEAD_WORM_STARVE_MASS,

    -- Count all worms died by old age
    (SELECT COUNT(*) FROM worms WHERE Notes LIKE "%old age%") AS DEAD_WORM_OLD_AGE_COUNT,

    -- Count mass of worms died by old age
    (SELECT SUM(Mass) FROM worms WHERE Notes LIKE "%old age%") AS DEAD_WORM_OLD_AGE_MASS,

    -- Count all dead mass
    (SELECT SUM(Mass) FROM worms WHERE Notes LIKE "%Cause of Death%") AS DEAD_MASS
