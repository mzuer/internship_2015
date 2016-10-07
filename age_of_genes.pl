#!/usr/bin/env perl -w


use lib "$ENV{HOME}/src/bioperl-1.6.1";
use lib "$ENV{HOME}/src/ensembl/modules";
use lib "$ENV{HOME}/src/ensembl-compara/modules";
use lib "$ENV{HOME}/src/ensembl-variation/modules";
use lib "$ENV{HOME}/src/ensembl-funcgen/modules";

use strict;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Registry;
$|=1;

# For the output, path to save the data
my $handle = undef;
my $fileMus = "/home/user/Bureau/ageMus.txt";  
my $fileXeno = "/home/user/Bureau/ageXeno.txt";
my $fileDanio = "/home/user/Bureau/ageDanio.txt";
my $fileCaeno = "/home/user/Bureau/ageCaeno.txt";
my $fileDroso = "/home/user/Bureau/ageDroso.txt";

# Connect to the Ensembl API
my $reg = "Bio::EnsEMBL::Registry";
$reg->load_registry_from_db(
  -host=>"ensembldb.ensembl.org",
  -user=>"anonymous",
);

# which Ensembl version is used?
use Bio::EnsEMBL::ApiVersion;
printf( "The API version used is %s\n", software_version() );

# Get all adaptors needed
my $gene_tree_adaptor = $reg->get_adaptor("Multi", "compara", "GeneTree");

# In a release database, there is a basal root that holds together all the trees.
my @children = @{$gene_tree_adaptor->fetch_all(-tree_type     => 'tree',
                                               -member_type   => 'protein',
                                               -clusterset_id => 'default',
                                           )}; # Protein tree only !

#my $counter=1; #counter for make a try before launch the all program

foreach my $tree (sort {$a->stable_id cmp $b->stable_id} (@children)) { #Sort to always get the same order

  # Tip: will make the traversal of the tree faster
  $tree->preload();

  # Stable ID of the tree (ENSGT...)
  print "Tree ", $tree->stable_id, "\n";
  # $tree->print_tree(10);

  # Get taxonomical level of the root
  my $tax_level = $tree->root->taxonomy_level(); 
  print "\tTaxonomy level: $tax_level\n";

  # Browse the tree and print for each species all the genes for each node

  foreach my $leaf (@{$tree->get_all_leaves}) {

    my $gene = $leaf->gene_member->stable_id;
    if ($gene =~m/ENSDARG/) { ## select zebrafish genes only
      # print on the screen      
      print $gene ,"\t", $tax_level, "\n";
      # print in a file
      open($handle, ">>", $fileDanio);
      print $handle $gene ,"\t", $tax_level, "\n";
      close $handle;
    } elsif($gene=~m/FBgn/) {## select drosophila genes only
      print $gene ,"\t", $tax_level, "\n";
      open($handle, ">>", $fileDroso); 
      print $handle $gene ,"\t", $tax_level, "\n";
      close $handle;      		
    } elsif($gene=~m/WBGene/){## select caenorhabditis genes only
      print $gene ,"\t", $tax_level, "\n";
      open($handle, ">>", $fileCaeno); 
      print $handle $gene ,"\t", $tax_level, "\n";
      close $handle;      
    } elsif($gene=~m/ENSXETG/){## select xenopus genes only
      print $gene ,"\t", $tax_level, "\n";
      open($handle, ">>", $fileXeno);
      print $handle $gene ,"\t", $tax_level, "\n";
      close $handle;      
    } elsif($gene=~m/ENSMUSG/){## select mouse genes only
      print $gene ,"\t", $tax_level, "\n";
      open($handle, ">>", $fileMus); 
      print $handle $gene ,"\t", $tax_level, "\n";
      close $handle;      
    }
  }
  $tree->root->release_tree;

#   last if ($counter++==10);   #counter for make a try before launch the all program

}
print "Done\n";
exit;

