while (<>) {
  ($east, $north, $elev, $mass) = split " ", $_;
  printf "$east $north $mass\n";
}
