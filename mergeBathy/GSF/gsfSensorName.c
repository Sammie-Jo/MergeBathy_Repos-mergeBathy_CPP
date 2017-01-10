#include "gsf.h"

char *gsfSensorName (gsfRecords *gsf)
{
   char *p = "unknown";
   switch (gsf->mb_ping.sensor_id)
   {
      case GSF_SWATH_BATHY_SUBRECORD_SEABEAM_SPECIFIC:
           p = "SeaBeam";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM12_SPECIFIC:
           p = "Simrad EM12";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM100_SPECIFIC:
           p = "Simrad EM100";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM950_SPECIFIC:
           p = "Simrad EM950";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM1000_SPECIFIC:
           p = "Simrad EM1000";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM121A_SPECIFIC:
           p = "Simrad EM121A";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM121_SPECIFIC:
           p = "Simrad EM121";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_SEAMAP_SPECIFIC:
           p = "SeaMap";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_SEABAT_SPECIFIC:
           if (gsf->mb_ping.sensor_data.gsfSeaBatSpecific.mode & 
               GSF_SEABAT_9002)
           {
              p = "Reson SeaBat 9002";
           }
           else if (gsf->mb_ping.sensor_data.gsfSeaBatSpecific.mode & 
               GSF_SEABAT_9003)
           {
              p = "Reson SeaBat 9003";
           }
           else
           {
              p = "Reson SeaBat 9001";
           }
           break;
      case GSF_SWATH_BATHY_SUBRECORD_SB_AMP_SPECIFIC:
           p = "SeaBeam (w/amp)";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_SEABAT_II_SPECIFIC:
           if (gsf->mb_ping.sensor_data.gsfSeaBatIISpecific.mode & 
               GSF_SEABAT_9002)
           {
              p = "Reson SeaBat 9002";
           }
           else if (gsf->mb_ping.sensor_data.gsfSeaBatSpecific.mode & 
               GSF_SEABAT_9003)
           {
              p = "Reson SeaBat 9003";
           }
           else
           {
              p = "Reson SeaBat 9001";
           }
           break;
      case GSF_SWATH_BATHY_SUBRECORD_SEABAT_8101_SPECIFIC:
           p = "Reson SeaBat 8101";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_RESON_8111_SPECIFIC:
           p = "Reson 8111";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_RESON_8124_SPECIFIC:
           p = "Reson 8124";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_RESON_8125_SPECIFIC:
           p = "Reson 8125";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_RESON_8150_SPECIFIC:
           p = "Reson 8150";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_RESON_8160_SPECIFIC:
           p = "Reson 8160";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_SEABEAM_2112_SPECIFIC:
           p = "SeaBeam 2112/36";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_ELAC_MKII_SPECIFIC:
           p = "ELAC MK II";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM300_SPECIFIC:
           p = "SIMRAD EM 300";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM1002_SPECIFIC:
           p = "SIMRAD EM 1002";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM3000_SPECIFIC:
           p = "SIMRAD EM 3000";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM120_SPECIFIC:
           p = "SIMRAD EM 120";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM122_SPECIFIC:
           p = "SIMRAD EM 122";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM710_SPECIFIC:
           p = "SIMRAD EM 710";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM3002_SPECIFIC:
           p = "SIMRAD EM 3002";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM3000D_SPECIFIC:
           p = "SIMRAD EM 3000D";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_EM3002D_SPECIFIC:
           p = "SIMRAD EM 3002D";
           break;
      case GSF_SINGLE_BEAM_SUBRECORD_ECHOTRAC_SPECIFIC:
      case GSF_SWATH_BATHY_SB_SUBRECORD_ECHOTRAC_SPECIFIC:
           p = "Echotrac";
           break;
      case GSF_SINGLE_BEAM_SUBRECORD_BATHY2000_SPECIFIC:
      case GSF_SWATH_BATHY_SB_SUBRECORD_BATHY2000_SPECIFIC:
           p = "Bathy 2000";
           break;
      case GSF_SINGLE_BEAM_SUBRECORD_MGD77_SPECIFIC:
      case GSF_SWATH_BATHY_SB_SUBRECORD_MGD77_SPECIFIC:
           p = "MGD77";
           break;
      case GSF_SINGLE_BEAM_SUBRECORD_BDB_SPECIFIC:
      case GSF_SWATH_BATHY_SB_SUBRECORD_BDB_SPECIFIC:
           p = "BDB";
           break; 
      case GSF_SINGLE_BEAM_SUBRECORD_NOSHDB_SPECIFIC:
      case GSF_SWATH_BATHY_SB_SUBRECORD_NOSHDB_SPECIFIC:
           p = "NOSHDB";
           break;
      case GSF_SWATH_BATHY_SB_SUBRECORD_PDD_SPECIFIC:
           p = "PDD";
           break;
      case GSF_SWATH_BATHY_SB_SUBRECORD_NAVISOUND_SPECIFIC:
           p = "NAVISOUND";
           break;
      case GSF_SWATH_BATHY_SUBRECORD_CMP_SASS_SPECIFIC:
           p = "Compressed SASS";
           break;
      default:
           p = "unknown";
   }
   return p;
}
