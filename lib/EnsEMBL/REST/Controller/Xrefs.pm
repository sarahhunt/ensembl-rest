package EnsEMBL::REST::Controller::Xrefs;
use Moose;
use namespace::autoclean;
use feature "switch";
require EnsEMBL::REST;
EnsEMBL::REST->turn_on_jsonp(__PACKAGE__);
use Try::Tiny;

BEGIN {extends 'Catalyst::Controller::REST'; }


# Expect 3 URLs coming from this class
#   - /xrefs/symbol/:species/:symbol
#   - /xrefs/id/:id
#   - /xrefs/name/:species/:name
# 
# Respond to params
#   - external_db='HGNC'
#   - db_type=core (might need to look in another DB)
#
# For /id/
#   - db_type=core
#   - object=gene
#   - species=human
#   - all_levels=1
#
# For /name

sub id_GET {}

sub id :Chained('/') PathPart('xrefs/id') Args(1)  ActionClass('REST') {
  my ($self, $c, $id) = @_;
  $c->stash()->{id} = $id;
  try {
    $c->log()->debug('Finding the object');
    my $obj = $c->model('Lookup')->find_object_by_stable_id($c, $id);
    $c->log()->debug('Processing the Xrefs');
    my $method = $c->request()->param('all_levels') ? 'get_all_DBLinks' : 'get_all_DBEntries';
    my $can = $obj->can($method);
    if(!$can) {
      my $msg = sprintf('The object type "%s" for ID "%s" cannot respond to the given request for Xrefs. Are you sure it has them?', ref($obj), $id);
      $c->log()->debug($msg);
      $c->go('ReturnError', 'custom', [$msg]);
    }
    my @args = ($obj);
    push(@args, $c->request->param('external_db')) if $c->request->param('external_db');
    my $entries = $can->(@args);
    $c->stash(entries => $entries);
    $c->forward('_encode');
  } 
  catch {
    $c->go('ReturnError', 'from_ensembl', [$_]);
  };
  $self->status_ok($c, entity => $c->stash->{entity});
}

sub symbol_GET {}

sub symbol :Chained('/') PathPart('xrefs/symbol') Args(2) ActionClass('REST') {
  my ($self, $c, $species, $symbol) = @_;
  $c->stash(species => $species, symbol => $symbol);
  my $external_db = $c->request->param('external_db');
  my $db_type = $c->request->param('db_type') || 'core';
  my @entries;
  try {
    my @objects_with_xrefs = $c->request->param('object') ? ($c->request->param('object')) : qw(gene transcript translation);
    foreach my $object_type (@objects_with_xrefs) {
      my $object_adaptor = $c->model('Registry')->get_adaptor($c->stash()->{species}, $db_type, $object_type);
      my $objects_linked_to_name = $object_adaptor->fetch_all_by_external_name($symbol, $external_db);
      while(my $obj = shift @{$objects_linked_to_name}) {
        my $encoded = {
          id => $obj->stable_id(),
          type => $object_type
        };
        push(@entries, $encoded);
      }
    }
  }
  catch {
    $c->go('ReturnError', 'from_ensembl', [$_]);
  };
  $self->status_ok( $c, entity => \@entries);
}

sub name_GET {}

sub name :Chained('/') PathPart('xrefs/name') Args(2) ActionClass('REST') {
  my ($self, $c, $species, $name) = @_;
  $c->stash(species => $species, name => $name);
  my $external_db = $c->request->param('external_db');
  my $db_type = $c->request->param('db_type') || 'core';
  try {
    my $dbentry_adaptor = $c->model('Registry')->get_adaptor($species, $db_type, 'dbentry');
    my $entries = $dbentry_adaptor->fetch_all_by_name($name, $external_db);
    $c->stash(entries => $entries);
    $c->forward('_encode');
  }
  catch {
    $c->go('ReturnError', 'from_ensembl', [$_]);
  };
  $self->status_ok($c, entity => $c->stash->{entity});
}

sub _encode :Private {
  my ($self, $c) = @_;
  my @encoded;
  foreach my $dbe (@{$c->stash()->{entries}}) {
    my $enc = {
      dbname          => $dbe->dbname(),
      db_display_name => $dbe->db_display_name(),
      display_id      => $dbe->display_id(),
      primary_id      => $dbe->primary_id(),
      description     => $dbe->description(),
      synonyms        => $dbe->get_all_synonyms(),
      version         => $dbe->version(),
      info_type       => $dbe->info_type(),
      info_text       => $dbe->info_text(),
    };
    given(ref($dbe)) {
      when('Bio::EnsEMBL::IdentityXref') {
        $enc->{xref_identity}     = $dbe->xref_identity();
        $enc->{xref_start}        = $dbe->xref_start();
        $enc->{xref_end}          = $dbe->xref_end();
        $enc->{ensembl_identity}  = $dbe->ensembl_identity();
        $enc->{ensembl_start}     = $dbe->ensembl_start();
        $enc->{ensembl_end}       = $dbe->ensembl_end();
        $enc->{score}             = $dbe->score();
        $enc->{evalue}            = $dbe->evalue();
        $enc->{cigar_line}        = $dbe->cigar_line() if $dbe->cigar_line();
      }
      when('Bio:EnsEMBL::OntologyXref') {
        $enc->{linkage_types} = $dbe->get_all_linkage_types();
      }
    }
    push(@encoded, $enc);
  }
  $c->stash(entity => \@encoded);
}

__PACKAGE__->meta->make_immutable;

1;
