#ifndef ASSERT
#define ASSERT(expr) if (!(expr)) {fprintf(stderr,"[ERROR](%s:%d) %s\n",__FILE__,__LINE__,#expr);exit(1);}
#endif
