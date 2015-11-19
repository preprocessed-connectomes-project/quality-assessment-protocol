require 'fileutils'
require 'rake'
require 'tmpdir'

destination = './docs/_site/'

desc 'Builds the website using Jekyll'
task :build do
   puts 'Rendering markdown to HTML with Jekyll'
   sh "cd docs && jekyll build"
   puts 'Rendering finished- results are in _site'
end

# From Evan Osenko @ https://evansosenko.com/posts/automatic-publishing-github-pages-travis-ci/
desc 'Generate deck from Travis CI and publish to GitHub Pages.'
task :travis do
  # if this is a pull request, do a simple build of the site and stop
  if ENV['TRAVIS_PULL_REQUEST'].to_s.to_i > 0
    puts 'Pull request detected. Executing build only.'
    sh 'bundle exec rake build'
    next
  end

  repo = %x(git config remote.origin.url).gsub(/^git:/, 'https:').strip
  deploy_url = repo.gsub %r{https://}, "https://#{ENV['GH_TOKEN']}@"
  deploy_branch = repo.match(/github\.io\.git$/) ? 'master' : 'gh-pages'
  rev = %x(git rev-parse HEAD).strip

  Dir.mktmpdir do |dir|
    dir = File.join dir, 'site'
    sh 'bundle exec rake build'
    fail "Build failed." unless Dir.exists? destination
    sh "git clone --branch #{deploy_branch} #{repo} #{dir}"
    sh %Q(rsync -rt --del --exclude=".git" --exclude=".nojekyll" #{destination} #{dir})
    Dir.chdir dir do
      # setup credentials so Travis CI can push to GitHub
      verbose false do
        sh "git config user.name '#{ENV['GIT_NAME']}'"
        sh "git config user.email '#{ENV['GIT_EMAIL']}'"
      end

      sh 'git add --all'
      sh "git commit -m 'Built from #{rev}'."
      verbose false do
        sh "git push -q #{deploy_url} #{deploy_branch}"
      end
    end
  end
end
