import pandas as pd


class CategoryGeneAnalyzer:
    def __init__(self, file_path):
        """
        Initialize the analyzer by loading category genes from a CSV file.

        Args:
        file_path (str): Path to the CSV file with UniProt and HUGO columns
        """
        # Read the CSV file
        self.gene_df = pd.read_csv(file_path)

        # Create a set of HUGO symbols for easy checking
        self.category_genes = set(self.gene_df['HUGO'])


    def validate_and_get_ratio(self, gene_set):
        """
        Validate if genes are in category genes and calculate ratio.

        Args:
        gene_set (set or list): Set or list of gene symbols to check

        Returns:
        tuple: (valid_genes, invalid_genes, ratio)
            - valid_genes: genes that are in category gene set
            - invalid_genes: genes that are not in category gene set
            - ratio: proportion of input genes that are category genes
        """
        # Convert input to a set to remove duplicates
        gene_set = set(gene_set)

        # Find valid and invalid genes
        valid_genes = gene_set.intersection(self.category_genes)
        invalid_genes = gene_set - self.category_genes

        # Calculate ratio
        ratio = len(valid_genes) / len(self.categorys_genes)

        return list(valid_genes), list(invalid_genes), ratio


# How to use it
# Path to the genes CSV file
file_path = ''


# Create an analyzer instance
analyzer = CategoryGeneAnalyzer(file_path)        

# Validation and ratio calculation
print("\nValidation and Ratio Check")
# Example gene sets to validate
test_sets = [
        ['PGK1', 'HK2', 'INVALID_GENE'],
        ['ENO1', 'PGK1', 'HK2']

for gene_set in test_sets:
  valid_genes, invalid_genes, ratio = analyzer.validate_and_get_ratio(gene_set)

  print(f"\nInput genes: {len(gene_set)}")
  print(f"Valid  genes: {len(valid_genes)}")
  print(f"Invalid genes: {len(invalid_genes)}")
  #print(f"Ratio of  genes: {ratio:.2f}")
  print(f"Valid genes recovered: {ratio*100:.2f}%")        