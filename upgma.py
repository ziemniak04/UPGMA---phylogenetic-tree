import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Union, Optional
from dataclasses import dataclass


@dataclass
class TreeNode:
    """Represents a node in the phylogenetic tree"""
    name: str
    left: Optional['TreeNode'] = None
    right: Optional['TreeNode'] = None
    distance: float = 0.0
    height: float = 0.0

    def is_leaf(self) -> bool:
        return self.left is None and self.right is None


class UPGMATreeBuilder:
    """
    UPGMA (Unweighted Pair Group Method with Arithmetic Mean) tree builder

    This class implements the UPGMA algorithm for constructing phylogenetic trees
    from either sequence data or distance matrices.
    """

    def __init__(self):
        self.sequences = {}
        self.distance_matrix = None
        self.sequence_names = []
        self.tree_root = None

    def load_sequences_from_fasta(self, fasta_content: str) -> bool:
        """
        Load sequences from FASTA format string

        Args:
            fasta_content: String containing FASTA formatted sequences

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            sequences = {}
            current_name = None
            current_seq = ""

            for line in fasta_content.strip().split('\n'):
                line = line.strip()
                if line.startswith('>'):
                    if current_name:
                        sequences[current_name] = current_seq.upper()
                    current_name = line[1:].strip()
                    current_seq = ""
                elif line:
                    current_seq += line

            if current_name:
                sequences[current_name] = current_seq.upper()

            if len(sequences) < 2:
                raise ValueError("At least 2 sequences required")

            lengths = [len(seq) for seq in sequences.values()]
            if len(set(lengths)) > 1:
                raise ValueError("All sequences must have the same length (MSA required)")

            self.sequences = sequences
            self.sequence_names = list(sequences.keys())
            return True

        except Exception as e:
            print(f"Error loading sequences: {e}")
            return False

    def load_sequences_from_list(self, sequences: List[Tuple[str, str]]) -> bool:
        """
        Load sequences from list of (name, sequence) tuples

        Args:
            sequences: List of (name, sequence) tuples

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            if len(sequences) < 2:
                raise ValueError("At least 2 sequences required")

            seq_dict = {}
            for name, seq in sequences:
                seq_dict[name] = seq.upper()

            lengths = [len(seq) for seq in seq_dict.values()]
            if len(set(lengths)) > 1:
                raise ValueError("All sequences must have the same length (MSA required)")

            self.sequences = seq_dict
            self.sequence_names = list(seq_dict.keys())
            return True

        except Exception as e:
            print(f"Error loading sequences: {e}")
            return False

    def load_distance_matrix(self, matrix: Union[np.ndarray, List[List[float]]],
                             names: List[str]) -> bool:
        """
        Load pre-computed distance matrix

        Args:
            matrix: Distance matrix as numpy array or list of lists
            names: List of sequence/species names

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            matrix = np.array(matrix)

            if matrix.shape[0] != matrix.shape[1]:
                raise ValueError("Distance matrix must be square")

            if len(names) != matrix.shape[0]:
                raise ValueError("Number of names must match matrix dimensions")

            if not np.allclose(matrix, matrix.T):
                raise ValueError("Distance matrix must be symmetric")

            if not np.allclose(np.diag(matrix), 0):
                raise ValueError("Diagonal elements must be zero")

            if np.any(matrix < 0):
                raise ValueError("Distance matrix cannot contain negative values")

            self.distance_matrix = matrix
            self.sequence_names = names.copy()
            return True

        except Exception as e:
            print(f"Error loading distance matrix: {e}")
            return False

    def calculate_distance_matrix(self) -> np.ndarray:
        """
        Calculate distance matrix from aligned sequences using Hamming distance

        Returns:
            np.ndarray: Symmetric distance matrix
        """
        if not self.sequences:
            raise ValueError("No sequences loaded")

        n = len(self.sequences)
        matrix = np.zeros((n, n))
        seq_list = [self.sequences[name] for name in self.sequence_names]

        for i in range(n):
            for j in range(i + 1, n):
                seq1, seq2 = seq_list[i], seq_list[j]
                differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')
                valid_positions = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')

                if valid_positions > 0:
                    distance = differences / valid_positions
                else:
                    distance = 1.0

                matrix[i][j] = matrix[j][i] = distance

        return matrix

    def build_upgma_tree(self) -> TreeNode:
        """
        Build phylogenetic tree using UPGMA algorithm

        Returns:
            TreeNode: Root node of the constructed tree
        """
        # Get or calculate distance matrix
        if self.distance_matrix is None:
            if not self.sequences:
                raise ValueError("No input data provided")
            distance_matrix = self.calculate_distance_matrix()
        else:
            distance_matrix = self.distance_matrix.copy()

        clusters = [TreeNode(name) for name in self.sequence_names]
        cluster_sizes = [1] * len(clusters)

        while len(clusters) > 1:
            n = len(clusters)

            min_dist = float('inf')
            min_i, min_j = -1, -1

            for i in range(n):
                for j in range(i + 1, n):
                    if distance_matrix[i][j] < min_dist:
                        min_dist = distance_matrix[i][j]
                        min_i, min_j = i, j

            new_node = TreeNode(
                name=f"({clusters[min_i].name},{clusters[min_j].name})",
                left=clusters[min_i],
                right=clusters[min_j],
                height=min_dist / 2
            )

            clusters[min_i].distance = new_node.height - clusters[min_i].height
            clusters[min_j].distance = new_node.height - clusters[min_j].height

            new_distances = []
            new_size = cluster_sizes[min_i] + cluster_sizes[min_j]

            for k in range(n):
                if k != min_i and k != min_j:
                    new_dist = (cluster_sizes[min_i] * distance_matrix[min_i][k] +
                                cluster_sizes[min_j] * distance_matrix[min_j][k]) / new_size
                    new_distances.append(new_dist)

            if min_j > min_i:
                clusters.pop(min_j)
                cluster_sizes.pop(min_j)
                clusters.pop(min_i)
                cluster_sizes.pop(min_i)
            else:
                clusters.pop(min_i)
                cluster_sizes.pop(min_i)
                clusters.pop(min_j)
                cluster_sizes.pop(min_j)

            clusters.append(new_node)
            cluster_sizes.append(new_size)

            new_matrix = np.zeros((len(clusters), len(clusters)))

            old_indices = [i for i in range(n) if i != min_i and i != min_j]
            for i, old_i in enumerate(old_indices):
                for j, old_j in enumerate(old_indices):
                    if i < j:
                        new_matrix[i][j] = new_matrix[j][i] = distance_matrix[old_i][old_j]

            for i, dist in enumerate(new_distances):
                new_matrix[i][-1] = new_matrix[-1][i] = dist

            distance_matrix = new_matrix

        self.tree_root = clusters[0]
        return self.tree_root

    def save_tree_newick(self, filename: str) -> bool:
        """
        Save tree in Newick format

        Args:
            filename: Output filename

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            if not self.tree_root:
                raise ValueError("No tree built yet")

            def node_to_newick(node: TreeNode) -> str:
                if node.is_leaf():
                    return f"{node.name}:{node.distance:.6f}"
                else:
                    left_str = node_to_newick(node.left)
                    right_str = node_to_newick(node.right)
                    return f"({left_str},{right_str}):{node.distance:.6f}"

            newick_str = node_to_newick(self.tree_root) + ";"

            with open(filename, 'w') as f:
                f.write(newick_str)

            return True

        except Exception as e:
            print(f"Error saving tree: {e}")
            return False

    def save_tree_table(self, filename: str) -> bool:
        """
        Save tree in tabular format

        Args:
            filename: Output filename

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            if not self.tree_root:
                raise ValueError("No tree built yet")

            with open(filename, 'w') as f:
                f.write("Node\tParent\tDistance\tHeight\tType\n")

                def write_node(node: TreeNode, parent_name: str = "ROOT"):
                    node_type = "Leaf" if node.is_leaf() else "Internal"
                    f.write(f"{node.name}\t{parent_name}\t{node.distance:.6f}\t{node.height:.6f}\t{node_type}\n")

                    if not node.is_leaf():
                        write_node(node.left, node.name)
                        write_node(node.right, node.name)

                write_node(self.tree_root)

            return True

        except Exception as e:
            print(f"Error saving tree table: {e}")
            return False

    def plot_tree(self, filename: str = None, show: bool = True) -> bool:
        """
        Create a clean horizontal dendogram visualization similar to standard phylogenetic tree plots

        Args:
            filename: Optional filename to save plot
            show: Whether to display the plot

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            if not self.tree_root:
                raise ValueError("No tree built yet")

            fig, ax = plt.subplots(1, 1, figsize=(10, 8))

            def get_leaves_ordered(node: TreeNode) -> List[str]:
                if node.is_leaf():
                    return [node.name]
                else:
                    return get_leaves_ordered(node.left) + get_leaves_ordered(node.right)

            leaves = get_leaves_ordered(self.tree_root)
            n_leaves = len(leaves)

            leaf_positions = {leaf: n_leaves - 1 - i for i, leaf in enumerate(leaves)}

            max_height = self.tree_root.height

            def draw_node(node: TreeNode) -> float:
                """
                Draw a node and return its y-position
                """
                if node.is_leaf():
                    y_pos = leaf_positions[node.name]

                    ax.plot([0, node.height], [y_pos, y_pos], 'k-', linewidth=2)

                    ax.text(-0.02 * max_height, y_pos, node.name,
                            ha='right', va='center', fontsize=11, fontweight='bold')

                    return y_pos
                else:
                    left_y = draw_node(node.left)
                    right_y = draw_node(node.right)

                    internal_y = (left_y + right_y) / 2.0

                    ax.plot([node.height, node.height], [left_y, right_y], 'k-', linewidth=2)

                    if node.left.height < node.height:
                        ax.plot([node.left.height, node.height], [left_y, left_y], 'k-', linewidth=2)
                    if node.right.height < node.height:
                        ax.plot([node.right.height, node.height], [right_y, right_y], 'k-', linewidth=2)

                    return internal_y

            draw_node(self.tree_root)

            margin_x_left = 0.05 * max_height
            margin_x_right = 0.05 * max_height
            margin_y = 0.3

            ax.set_xlim(-margin_x_left, max_height + margin_x_right)
            ax.set_ylim(-margin_y, n_leaves - 1 + margin_y)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_yticks([])

            ax.set_xlabel('Genetic Distance', fontsize=13, fontweight='bold')
            ax.spines['bottom'].set_linewidth(1.5)

            n_ticks = 6
            x_ticks = np.linspace(0, max_height, n_ticks)
            ax.set_xticks(x_ticks)
            ax.set_xticklabels([f'{x:.3f}' for x in x_ticks], fontsize=10)

            ax.grid(True, axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
            ax.set_axisbelow(True)

            ax.set_title('UPGMA Phylogenetic Tree (Dendogram)',
                         fontsize=16, fontweight='bold', pad=20)

            ax.set_facecolor('white')
            fig.patch.set_facecolor('white')

            plt.tight_layout()

            if filename:
                plt.savefig(filename, dpi=300, bbox_inches='tight',
                            facecolor='white', edgecolor='none')
                print(f"Tree plot saved to {filename}")

            if show:
                plt.show()
            else:
                plt.close()

            return True

        except Exception as e:
            print(f"Error plotting tree: {e}")
            return False


def get_user_input_sequences():
    """
    Get sequence data from user input
    Returns:
        List[Tuple[str, str]]: List of (name, sequence) tuples
    """
    sequences = []
    print("\nEnter aligned sequences (MSA format):")
    print("Type 'done' when finished entering sequences")
    print("Format: name sequence (separated by space)")
    print("Example: Human ATGCGTACGTACGT")

    while True:
        user_input = input("\nEnter sequence (or 'done'): ").strip()

        if user_input.lower() == 'done':
            break

        try:
            parts = user_input.split(None, 1)
            if len(parts) != 2:
                print("Error: Please enter name and sequence separated by space")
                continue

            name, sequence = parts
            sequence = sequence.replace(' ', '').upper()

            if not sequence:
                print("Error: Empty sequence")
                continue

            if not all(c in 'ATGCN-' for c in sequence):
                print("Warning: Sequence contains non-standard nucleotide characters")

            sequences.append((name, sequence))
            print(f"✓ Added: {name} ({len(sequence)} bp)")

        except Exception as e:
            print(f"Error parsing input: {e}")

    return sequences


def get_user_input_distance_matrix():
    """
    Get distance matrix data from user input

    Returns:
        Tuple[List[List[float]], List[str]]: Distance matrix and names
    """
    print("\nEnter distance matrix:")

    # Get number of species
    while True:
        try:
            n = int(input("Number of species/taxa: "))
            if n < 2:
                print("Error: At least 2 species required")
                continue
            break
        except ValueError:
            print("Error: Please enter a valid number")

    names = []
    print(f"\nEnter {n} species/taxa names:")
    for i in range(n):
        while True:
            name = input(f"Species {i + 1} name: ").strip()
            if name:
                names.append(name)
                break
            print("Error: Name cannot be empty")

    print(f"\nEnter distance matrix ({n}x{n}):")
    print("Enter each row on a separate line, with values separated by spaces")
    print("Note: Matrix should be symmetric with zeros on diagonal")

    matrix = []
    for i in range(n):
        while True:
            try:
                row_input = input(f"Row {i + 1} ({names[i]}): ").strip()
                row = [float(x) for x in row_input.split()]

                if len(row) != n:
                    print(f"Error: Expected {n} values, got {len(row)}")
                    continue

                matrix.append(row)
                break

            except ValueError:
                print("Error: Please enter numeric values separated by spaces")

    return matrix, names


def main():
    """
    Main function with interactive input for UPGMA tree builder
    """
    print("UPGMA Phylogenetic Tree Constructor")
    print("=" * 40)
    print("\nChoose input method:")
    print("1. Enter aligned sequences (MSA)")
    print("2. Enter distance matrix")
    print("3. Run examples")

    while True:
        choice = input("\nEnter your choice (1-3): ").strip()

        if choice == "1":
            # Sequence input mode
            print("\n" + "=" * 50)
            print("SEQUENCE INPUT MODE")
            print("=" * 50)

            sequences = get_user_input_sequences()

            if len(sequences) < 2:
                print("Error: At least 2 sequences required")
                return

            builder = UPGMATreeBuilder()

            if builder.load_sequences_from_list(sequences):
                print("\n✓ Sequences loaded successfully")

                tree = builder.build_upgma_tree()
                print("✓ UPGMA tree constructed")

                prefix = input("\nEnter output filename prefix (default: 'tree'): ").strip()
                if not prefix:
                    prefix = "tree"

                builder.save_tree_newick(f"{prefix}.newick")
                builder.save_tree_table(f"{prefix}.txt")
                print(f"✓ Tree saved to {prefix}.newick and {prefix}.txt")

                plot_choice = input("Generate tree plot? (y/n): ").strip().lower()
                if plot_choice in ['y', 'yes']:
                    builder.plot_tree(f"{prefix}.png", show=True)
                    print(f"✓ Tree plot saved to {prefix}.png")

            break

        elif choice == "2":
            print("\n" + "=" * 50)
            print("DISTANCE MATRIX INPUT MODE")
            print("=" * 50)

            matrix, names = get_user_input_distance_matrix()

            builder = UPGMATreeBuilder()

            if builder.load_distance_matrix(matrix, names):
                print("\n✓ Distance matrix loaded successfully")

                tree = builder.build_upgma_tree()
                print("✓ UPGMA tree constructed")

                prefix = input("\nEnter output filename prefix (default: 'tree'): ").strip()
                if not prefix:
                    prefix = "tree"

                builder.save_tree_newick(f"{prefix}.newick")
                builder.save_tree_table(f"{prefix}.txt")
                print(f"✓ Tree saved to {prefix}.newick and {prefix}.txt")

                plot_choice = input("Generate tree plot? (y/n): ").strip().lower()
                if plot_choice in ['y', 'yes']:
                    builder.plot_tree(f"{prefix}.png", show=True)
                    print(f"✓ Tree plot saved to {prefix}.png")

            break

        elif choice == "3":
            print("\n" + "=" * 50)
            print("RUNNING EXAMPLES")
            print("=" * 50)

            print("\nExample 1: Building tree from sequence data")
            print("-" * 45)

            sample_sequences = [
                ("Human", "ATGCGTACGTACGT"),
                ("Chimp", "ATGCGTACGTACGA"),
                ("Mouse", "ATGCGTACGTACCG"),
                ("Rat", "ATGCGTACGTACCA")
            ]

            builder1 = UPGMATreeBuilder()

            if builder1.load_sequences_from_list(sample_sequences):
                print("✓ Sequences loaded successfully")

                tree = builder1.build_upgma_tree()
                print("✓ UPGMA tree constructed")

                builder1.save_tree_newick("example1_tree.newick")
                builder1.save_tree_table("example1_tree.txt")
                print("✓ Tree saved to files")

                builder1.plot_tree("example1_tree.png", show=False)
                print("✓ Tree plot saved")

            print("\nExample 2: Building tree from distance matrix")
            print("-" * 48)

            distance_matrix = [
                [0.0, 0.2, 0.6, 0.8],
                [0.2, 0.0, 0.7, 0.9],
                [0.6, 0.7, 0.0, 0.3],
                [0.8, 0.9, 0.3, 0.0]
            ]

            species_names = ["Species_A", "Species_B", "Species_C", "Species_D"]

            builder2 = UPGMATreeBuilder()

            if builder2.load_distance_matrix(distance_matrix, species_names):
                print("✓ Distance matrix loaded successfully")
                tree = builder2.build_upgma_tree()
                print("✓ UPGMA tree constructed")

                builder2.save_tree_newick("example2_tree.newick")
                builder2.save_tree_table("example2_tree.txt")
                print("✓ Tree saved to files")

                builder2.plot_tree("example2_tree.png", show=False)
                print("✓ Tree plot saved")

            print("\nExample analysis complete! Check the generated files:")
            print("- example1_tree.newick, example1_tree.txt, example1_tree.png")
            print("- example2_tree.newick, example2_tree.txt, example2_tree.png")

            break

        else:
            print("Invalid choice. Please enter 1, 2, or 3.")


if __name__ == "__main__":
    main()