
# Note that need to hit "enter" twice to set blank password during mysql-server install

apt-get install -y libmysqlclient-dev mysql-server
cpanm cpanm DBD::mysql
