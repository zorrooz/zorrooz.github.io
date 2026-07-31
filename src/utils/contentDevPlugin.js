/**
 * Vite dev plugin: regenerate content when content-src changes.
 * New/removed files require a server restart (import.meta.glob is resolved
 * at module transform time and would otherwise miss them).
 */

import { watch } from 'chokidar'

import { runAllGenerators } from './runAllGenerators.js'

export function contentDev() {
  let timer = null
  return {
    name: 'gblog:content-dev',
    apply: 'serve',
    async configureServer(server) {
      const run = async (label) => {
        try {
          await runAllGenerators()
          server.config.logger.info(`[content] regenerated (${label})`)
          server.ws.send({ type: 'full-reload' })
        } catch (err) {
          server.config.logger.error(`[content] regeneration failed: ${err.message}`)
        }
      }

      await run('startup')

      const watcher = watch('src/content-src', { ignoreInitial: true })
      watcher.on('all', (event) => {
        clearTimeout(timer)
        timer = setTimeout(() => {
          const needsRestart = event === 'add' || event === 'unlink'
          run(needsRestart ? 'restart' : 'change').then(() => {
            if (needsRestart) server.restart()
          })
        }, 300)
      })

      server.httpServer?.on('close', () => watcher.close())
    }
  }
}
