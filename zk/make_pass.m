%quick and dirty
%from cellsArray get a unit meta and run the visualize response
a_unit=cellsArray(3)

vr=visualize_responses_play(a_unit.mouse,a_unit.sess,a_unit.rec,a_unit.clu,'odor')
