struct cmp_str
{
  bool operator()(char const *a, char const *b) const
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef std::map<char*,std::vector<int>,cmp_str> mm;
mm* read_pops(char *fname);
