% from sources
standard_channels = { ...
	'satellites'								'count'; ...
	'time'										's'; ...
	'latitude'									'minutes'; ...
	'longitude'									'minutes'; ...
    'latitude_raw'                              'minutes'; ...
	'longitude_raw'								'minutes'; ...
	'velocity'									'knots'; ...
	'velocity knots'							'knots'; ...
	'velocity kmh'								'km/h'; ...
	'heading'									'deg'; ...
	'height'									'm'; ...
    'height_raw'                                'm'; ...
	'vertical velocity m/s'						'm/s'; ...
	'vertical velocity'							'm/s'; ...
	'velocity quality'							''; ...
	'long accel g'								'g'; ...
	'lat accel g'								'g'; ...
	'brake distance m'							'm'; ...
	'distance m'								'm'; ...
	'adc01-1'									''; ...
	'adc01-2'									''; ...
	'adc01-3'									''; ...
	'adc01-4'									''; ...
	'glonass_sats'								'count'; ...
	'gps_sats'									'count'; ...
	'yaw rate'									'deg/s'; ...
	'latacc'									'm/s^2'; ...
	'reserved'									''; ...
	'yaw rate 2'								'deg/s'; ...
	'latacc 2'									'm/s^2'; ...
	'solution type'								''; ...
	'velocity quality'							''; ...
	'internal temperature'						'deg'; ...
	'cf_buffersize'								''; ...
	'ramaddress'								''; ...
	'event 1 time'								's'; ...
	'event 2 time'								's'; ...
	'battery 1 voltage'							'V'; ...
	'battery 2 voltage'							'V'; ...
% from tex's list
	'absolute heading'							'deg'; ...
	'gps_longacc'								''; ...
	'gps_latacc'								''; ...
	'yaw rate deg/s'							'deg/s'; ...
	'yaw_rate2'									'deg/s'; ...
	'yaw_rate'									'deg/s'; ...
	'yawrate'									'deg/s'; ...
	'lateral accel (yaw sensor)'				''; ...
	'latacc_2'									''; ...
	'slip angle'								'deg'; ...
	'slip_angle'								'deg'; ...
	'radius of turn'							'm'; ...
	'trigger event time in clock counts'		''; ...
	'distance feet'								''; ...
	'gps distance from start point in metres'	'm'; ...
	'gps distance from start point in feet'		'ft'; ...
	'incremental time in seconds'				's'; ...
	'long acc calc'								'm/s^2'; ...
	'lat acc calc'								'm/s^2'; ...
	'x position in'								''; ...
	'y position in'								''; ...
	'relative height'							'm'; ...
	'relative height feet'						'ft'; ...
	'centre line deviation'						'm'; ...
    'roll_imu'                                  'deg'; ...
    'roll_angle'                                'deg/s'; ...
    'pitch_ang.'                                'deg'; ...
    'slip_fl'                                   'deg'; ...
    'slip_fr'                                   'deg'; ...
    'slip_rl'                                   'deg'; ...
    'slip_rr'                                   'deg'; ...
    'slip_cog'                                  'deg'; ...
    'x_accel'                                   'g'; ...
    'y_accel'                                   'g'; ...
    'z_accel'                                   'g'; ...
    'rollrate'                                  'deg/s'; ...
    'true_head2'                                'deg'; ...
    'head_imu'                                  'deg'; ...
    'heading_raw'                               'deg'; ...
    'pitch_imu'                                 'deg'; ...
    'pos.qual.'                                 'count'; ...
    'lng_jerk'                                  'm/s^^2'; ...
    'lat_jerk'                                  'm/s^^2'; ...
    'head_imu2'                                 'deg'; ...
    'vertical velocity kmh'                     'km/h'; ...
    'vertical_velocity_raw'                     'km/h'; ...
    'speed_raw'                                 'km/h'; ...
    % from S90 CAN bus
    'accrpedlln'                                '%'; ...
    'AccrPedlLnr'                               '%'; ...
    'engn'                                      'rpm'; ...
    'ptgearact'                                 'count'; ...
    'pttqatw_fl'                                'Nm'; ...
    'pttqatw_fr'                                'Nm'; ...
    'swa'                                       'rad'; ...
    'swaspd'                                    'rad/s'; ...
    'flwhlspd'                                  'm/s'; ...
    'frwhlspd'                                  'm/s'; ...
    'rlwhlspd'                                  'm/s'; ...
    'rrwhlspd'                                  'm/s'; ...
    'brkpedlpsd'                                '0 = off, 1 = on'; ...
    'vehspdlgt'                                 'm/s'; ...
    'algt1'                                     'm/s^2'; ...
    'alat1'                                     'm/s^2'; ...
    'rollrate1'                                 'rad/s'; ...
    'yawrate1'                                  'rad/s'; ...
    'temp'                                      'deg'; ...
};