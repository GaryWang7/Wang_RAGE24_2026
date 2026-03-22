script=<script_dir>

singularity exec \
--cleanenv \
--bind $HOME \
--bind /project/parkercwlab \
$HOME/images/isoseq_1.0.4.sif \
bash $script


