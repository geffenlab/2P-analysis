
% update suite2P_make_db_KCW with your mouse info
% run suite2P_master_file_KCW
% run new_main, sort neurons, and save
% processSingleBlocks OR separateBlocksSortedTogether if there is more than
% one expt sorted together (i.e. multiple folders in your expt folders)
% run makeBasicFile

%%FOR RED CHANNEL - MT procedure = I have a two channel recording of 1min
%%silence and green-channel recordings with stims and all
% update suite2P_make_db_KCW with the green channel recordings only
% run suite2P_master_file_KCW
% run new_main, sort neurons, and save
% processSingleBlocks OR separateBlocksSortedTogether if there is more than
% one expt sorted together (i.e. multiple folders in your expt folders)
% run makeBasicFile
%%%% 
% then for red cell identification
%  update suite2P_make_db_KCW with the two channel recording of silence
% run suite2P_master_file_KCW_MT_RedChannel: this 
% transfer the registered tiffs to your ocmputer
% concatenate in ImageJ the green channels tiffs together, save ; make a
% z-project average of the whole recording, save;
% do the same with the red channel tiffs
% run the program Indentification_RedCell_Movie_MT ... something

