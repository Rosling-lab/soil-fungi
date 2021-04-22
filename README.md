## Soil fungi bioinformatics pipeline
This repository hosts a generic pipeline to make ASVs from raw PacBio reads.

### To run on UPPMAX

* Check out the repository into a new project directory (replace `<proj_directory>` with a name for your project):
  ```bash
  git clone git@github.com:Rosling-lab/soil-fungi.git <proj_directory>
  ```
* Inside the new directory, link to the raw data (replace `<raw_data_path>` with the path to your raw data, e.g., `/proj/rosling_storage/Archy/b2013211/community/raw-data/pb_363`):
  ```bash
  cd <proj_directory>
  mkdir rawdata
  cd rawdata
  ln -s <raw_data_path> .
  ```
  (You can repeat the last line several times if you want to analyze several sequencing runs.)
* Schedule the script to run:
  ```bash
  sbatch run-node.sh
  ```

### Pipeline steps

- [X] Run improved CCS on RSII reads.
- [ ] Demultiplexing
- [X] Orient reads and trim primers (ITS1 and LR5).
  - [ ] Option for different primer pairs.
- [X] Filter: min length 1000 bp, max length 2000 bp, max 1% expected errors.
  - [ ] Option for different filter criteria.
- [ ] Generate consensus ASVs with LSUx + dada2 + Tzara
- [ ] Assign taxonomy with Silva/RDP/Unite
- [ ] Generate a tree
- [ ] Refine taxonomy with PHYLOTAX
