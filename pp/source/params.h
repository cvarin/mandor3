//main parameters for the entire program

#define MAXSTRLEN 1000		//maximum supposed string length in the cfg file
#define STRLEN 256		//supposed string length for file names/error messages/etc.
#define VERT_EXTRA 50		//extra pixels added to vert_pic_size to account for height of title(s)
#define CPU_64_BIT 1		//whether the machine is 64-bit (1) or a 32-bit (0) one
#define SMART_READING 0		//whether skip unnecessary data bit while reding binary files or not
#define MAX_Q_DIV_M 100		//if q/m in config file is larger than this value, process all the q/ms
#define VERY_SMALL_VALUE 1E-100	//very small value
#define SMALL_VALUE 1E-4	//small value
#define PI 3.141592653589793238462643
#define WEAKEN_COARSENING 1	//factor by which the coarsening of meshed data is weakened (i.e. WEAKEN_COARSENING*pic_size points are saved rather than just pic_size) 
#define BUFFER_SIZE 1048576	//size of the buffer for various intermediate tasks

//config files
#define CFG_COMMON "configs/common.cfg"
#define CFG_MDATA "configs/meshed_data.cfg"
#define CFG_MCUTS "configs/meshed_data.cfg"
#define CFG_PS "configs/scattered_data.cfg"
#define CFG_DF "configs/scattered_data.cfg"
#define CFG_ST "configs/special_tasks.cfg"


