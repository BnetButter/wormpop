SELECT
    (SELECT SUM(Egg_to_Larva) FROM stage_transition) as Egg_to_Larva_CNT,
    (SELECT SUM(Egg_to_Larva_Mass) FROM stage_transition) as Egg_to_Larva_Mass,
    (SELECT SUM(Dauer_to_Larva) FROM stage_transition) as Dauer_to_Larva_CNT,
    (SELECT SUM(Dauer_to_Larva_Mass) FROM stage_transition) as Dauer_to_Larva_Mass,
    (SELECT SUM(l1arrest_to_larva) FROM stage_transition) as l1arrest_to_larva_CNT,
    (SELECT SUM(l1arrest_to_larva_mass) FROM stage_transition) as l1arrest_to_larva_MASS,
    (SELECT SUM(Larva_culled_ind) FROM dynamic_stage_transition) as Larva_Culled_CNT,
    (SELECT SUM(Larva_culled_mass) FROM dynamic_stage_transition) as Larva_Culled_MASS,
    (SELECT SUM(Larva_starvation_ind) FROM dynamic_stage_transition) as Larva_starvation_CNT,
    (SELECT SUM(Larva_starvation_mass) FROM dynamic_stage_transition) as Larva_starvation_MASS,
    (SELECT SUM(Larva_arrested_development_ind) FROM dynamic_stage_transition) as Larva_starvation_AD_CNT,
    (SELECT SUM(Larva_arrested_development_mass) FROM dynamic_stage_transition) as Larva_starvation_MASS_CNT;