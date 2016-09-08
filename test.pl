#!/usr/bin/perl -w
use strict;
use warnings;
use Text::Mining;
use Text::Mining::Parser::PubMedCentral;

    my $tm = Text::Mining->new();

    my $corpus = $tm->get_corpus({ corpus_id => 1 });
    my $result = $corpus->add_dir({ dir    => '/home/user/data', 
                                    parser => 'PubMedCentral' });