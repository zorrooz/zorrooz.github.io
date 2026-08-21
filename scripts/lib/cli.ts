/**
 * CLI 引导共享实现（自 scripts/generators/core 按职责拆分迁移）。
 * 行为与迁移前逐字节一致；core/index.ts 仅 re-export，消费方 import 路径不变。
 */

import path from 'path'
import { fileURLToPath } from 'url'

/**
 * 判断模块是否被直接执行（node foo.ts），容忍相对/绝对 argv 形式。
 * 用于门控 CLI 入口，避免 import 生成器模块时误触发其 main()。
 */
export function isDirectRun(importMeta: ImportMeta): boolean {
  return Boolean(process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(importMeta.url))
}

/** 单个生成器 CLI 的通用引导：banner + 执行 + 失败退出码 */
export function runCliScript(name: string, run: () => void | Promise<void>): void {
  console.log(`Starting ${name} generation script...`)
  const done = () => console.log(`${name} generation complete.`)
  try {
    const result = run()
    if (result instanceof Promise) {
      result.then(done).catch((err) => {
        console.error(`${name} generation failed:`, err)
        process.exitCode = 1
      })
    } else {
      done()
    }
  } catch (error) {
    console.error(`${name} generation failed:`, error)
    process.exitCode = 1
  }
}

/** 统一的「生成成功」日志，可附括号明细 */
export function logWriteSuccess(targetPath: string, detail?: string): void {
  console.log(`Successfully generated: ${targetPath}${detail ? ` (${detail})` : ''}`)
}
