%for debugging trial number and trial pin lookup
mouse = 'ZKawakeM72';
sess  = 7;
rec   = 'c';
run  = 1;
ds=data_structure_tools;

event = struct('event',[],'type',[],'chanId',[],'evtFcn',[]);
eventsList = [];

% trial pin
event.event  = 'trialPin';
event.type   = 'digital';
event.chanId = 'trPin';
event.plotcolor    = 'r';
event.evtFcn = @get_analog_events;
eventsList   = [eventsList event];
 

% final valve pin
event.event  = 'finalValve';
event.type   = 'digital';
event.chanId = 'FVPin';
event.plotcolor     = 'g';
event.evtFcn = @get_analog_events;
eventsList   = [eventsList event];

% laser
event.event  = 'laser';
event.type   = 'semiDigital';
event.chanId = 'Laser';
event.plotcolor    = 'b';
event.evtFcn = @get_analog_events;
eventsList   = [eventsList event];


%get all the events (streams and events)

for ie=1:numel(eventsList)
    event=eventsList(ie);
    evName = event.event;
    evt.(evName).stream = ds.get_data_stream(event.chanId,mouse,sess,rec,run);
    evt.(evName).events = ds.get_analog_events(event.chanId,mouse,sess,rec,run,'figures','noplot');
    plot(evt.(evName).stream.data(pr),event.plotcolor);
end
evt.trialNum.stream = ds.get_data_stream('trNum',mouse,sess,rec,run);
evt.trialNum.events = get_trial_numbers(mouse,sess,rec,run);

figure(1)
hold on
pr = 1:2000000;
for ie=1:numel(eventsList)
    event=eventsList(ie);
    evName = event.event;
    plot(evt.(evName).stream.data(pr),event.plotcolor);
end
plot(32768-evt.trialNum.stream.data(pr),'k');

figure(2)
hold on
pr = numel(evt.(evName).stream.data) + (-2000000:0);
for ie=1:numel(eventsList)
    event=eventsList(ie);
    evName = event.event;
    plot(evt.(evName).stream.data(pr),event.plotcolor);
end
plot(32768-evt.trialNum.stream.data(pr),'k');

figure(3)
hold on
pr = (8363088):(11460528);
for ie=1:numel(eventsList)
    event=eventsList(ie);
    evName = event.event;
    plot(evt.(evName).stream.data(pr),event.plotcolor);
end
plot(32768-evt.trialNum.stream.data(pr),'k');