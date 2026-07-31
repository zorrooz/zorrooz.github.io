/**
 * Vite dev plugin: regenerate content when content-src changes.
 * New/removed files require a server restart (import.meta.glob is resolved
 * at module transform time and would otherwise miss them).
 */

import { watch } from 'chokidar'
import path from 'path'
import type { Plugin, ViteDevServer } from 'vite'

import { runAllGenerators } from './runAllGenerators.ts'

const contentSrcDir = path.resolve(import.meta.dirname, '../content-src')

export function contentDev(): Plugin {
  let timer: ReturnType<typeof setTimeout> | null = null
  return {
    name: 'gblog:content-dev',
    apply: 'serve',
    async configureServer(server: ViteDevServer) {
      const run = async (label: string) => {
        try {
          await runAllGenerators()
          server.config.logger.info(`[content] regenerated (${label})`)
          server.ws.send({ type: 'full-reload' })
        } catch (err) {
          server.config.logger.error(`[content] regeneration failed: ${(err as Error).message}`)
        }
      }

      await run('startup')

      const watcher = watch(contentSrcDir, { ignoreInitial: true })
      watcher.on('all', (event) => {
        if (timer) clearTimeout(timer)
        timer = setTimeout(() => {
          const needsRestart = event === 'add' || event === 'unlink'
          run(needsRestart ? 'restart' : 'change').then(() => {
            if (needsRestart) server.restart()
          })
        }, 300)
      })

      server.httpServer?.on('close', () => watcher.close())
    },
  }
}
