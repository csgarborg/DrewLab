eq='motVecY_rosPos ~ xPos+yPos+sex+ (1|mouseID)';%+weight 
glm_model=fitglme( glmTableData, eq);