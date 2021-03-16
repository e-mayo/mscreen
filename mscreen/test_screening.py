   class Test_Vina:

        ligands = Path("../../data/ligands")
        receptors = Path("../../data/receptor")
        out = Path("../../data/out")
        conf = Path("../../data/conf.txt")
        backend = "vina"

        def test_prepared_folder(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend)
            s.prepared_folder()

        def test_prepare_receptors(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend)
            s.prepared_folder()
            s.prepare_receptors()

        def prepare_ligands(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend)
            s.prepared_folder()
            s.prepare_ligands()

        def prepare_screening(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend)
            s.prepare_screening()

        def test_get_docking_executable(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend)
            docking_engine = {
                "fwavina",
                "vina",
                "gwovina",
                "ledock",
                "psovina",
                "qvina2",
                "qvina-w",
                "smina",
            }
            for backend in docking_engine:
                s.get_docking_executable(backend)

        def test_run_vina(self):
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend)

            s.run_vina(self, 'out_path', 'lig_path', 'rec_path')

        def test_run_screening_prepare_true():
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend, prepare=True)
            s.run_screening()

        def test_run_screening_prepare_False():
            s = VinaScreening(self.ligands, self.receptors,
                              self.out, self.backend, prepare=False)
            s.run_screening()
