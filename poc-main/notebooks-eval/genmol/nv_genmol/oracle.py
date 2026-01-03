class Oracle:

    @staticmethod
    def TdcScore(name):
        try:
            import tdc
        except:
            print('This oracle requires the pyTDC package. Please run "pip install pyTDC" to install it.')
            exit()
        
        return tdc.Oracle(name)

    @staticmethod
    def RDKitScore(name):
        from rdkit import Chem

        def ScoreFromSmiles(func):
            def score(smiles):
                try: 
                    mol_obj = Chem.MolFromSmiles(smiles)
                except:
                    mol_obj = None

                if mol_obj is None:
                    return None
                    
                return func(mol_obj)
            
            return score

        if name == 'QED':
            from rdkit.Chem.QED import qed
            return ScoreFromSmiles(qed)
        elif name == 'LogP':
            from rdkit.Chem.Crippen import MolLogP
            return ScoreFromSmiles(MolLogP)
        elif name == 'SA':
            from rdkit.Contrib.SA_Score import sascorer
            return ScoreFromSmiles(sascorer.calculateScore)
        else:
            assert False, "Unsupported score name: " + name
    
    @staticmethod
    def HSECScore(score_type='hsec', ref_smiles=None):
        """
        scorer using the MoleculeScorer for HSEC/AIH scoring.
        
        Parameters:
        score_type (str): Type of score to return - 'hsec', 'aih', or 'both'
        
        Returns:
        function: A scoring function that accepts a SMILES string and returns the requested score
        """
        # Import the MoleculeScorer (assuming it's from a module)
        from nv_genmol.by_scorer_new import MoleculeScorer  # Replace 'your_module' with actual module name
        import pandas as pd
        
        # Create the scorer with your specifications
        hsec_scorer = MoleculeScorer.create_scorer(
            reference_smiles=ref_smiles,
            use_qsar_hsec=True,
            use_qsar_aih=False,
            use_qed=False,
            use_similarity=True,
            similarity_weight=1.0,
            similarity_threshold=1.0,
        )
        
        # Define function that handles a single SMILES string
        def score_single_smiles(smiles):
            
            # Score the molecule
            result = hsec_scorer.score([smiles])
            
            # Return the appropriate score based on score_type
            if score_type == 'hsec':
                return result.qsar_hsec_scores[0]
            if score_type == 'comp':
                pos_res = -1.0 * result.combined_scores[0]
                return pos_res
                
        return score_single_smiles
    
    
    def __init__(self, score = None):
        self.score = score

    def evaluate(self, molecules):
        return {_:self.score(_) for _ in molecules}



