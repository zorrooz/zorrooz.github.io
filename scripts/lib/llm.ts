// scripts/lib/llm.ts — LLM 配置 / 客户端 / 调用骨架的唯一实现（translator 与 tagMerger 共用）
import OpenAI from 'openai'
import fs from 'fs'
import path from 'path'
import { pathToFileURL } from 'url'

/** llmConfig.ts（gitignore）支持的字段；thinking 缺省 = 关闭思考模式 */
export interface LlmConfig {
  url: string
  apikey: string
  model: string
  /** 是否启用 DeepSeek 思考模式（默认关闭：省 token 且恢复 temperature 语义） */
  thinking?: boolean
}

let clientCache: OpenAI | null = null

/**
 * 懒加载 DeepSeek 配置（llmConfig.ts 被 gitignore，CI 无此文件）。
 * 模块可能被 vite-ssg 打包，静态 import 会因缺文件失败，故动态 import。
 */
export async function loadLlmConfig(): Promise<LlmConfig | null> {
  const cfgPath = path.join(import.meta.dirname, '../llmConfig.ts')
  if (!fs.existsSync(cfgPath)) return null
  try {
    const mod = (await import(pathToFileURL(cfgPath).href)) as {
      default?: Partial<LlmConfig>
    }
    const cfg = mod.default || {}
    if (!cfg.url || !cfg.apikey || !cfg.model) return null
    return { url: cfg.url, apikey: cfg.apikey, model: cfg.model, thinking: cfg.thinking }
  } catch {
    return null
  }
}

export async function getLlmClient(): Promise<OpenAI> {
  if (clientCache) return clientCache
  const cfg = await loadLlmConfig()
  if (!cfg)
    throw new Error('未找到 LLM 配置：请在 scripts/llmConfig.ts 配置 { url, apikey, model[, thinking] }')
  clientCache = new OpenAI({ baseURL: cfg.url, apiKey: cfg.apikey })
  return clientCache
}

export interface CompleteChatOptions {
  temperature?: number
  /** 覆盖配置的 thinking 开关（缺省跟随 llmConfig.ts） */
  thinking?: boolean
}

export interface CompleteChatResult {
  text: string
  usage?: { promptTokens?: number; completionTokens?: number; totalTokens?: number }
}

/** 统一的 chat.completions 调用骨架：返回首条 assistant 消息文本与 token 用量 */
export async function completeChat(
  systemPrompt: string,
  userContent: string,
  options: CompleteChatOptions = {},
): Promise<CompleteChatResult> {
  const [client, cfg] = await Promise.all([getLlmClient(), loadLlmConfig()])
  const thinking = options.thinking ?? cfg?.thinking === true
  const completion = await client.chat.completions.create({
    messages: [
      { role: 'system', content: systemPrompt },
      { role: 'user', content: userContent },
    ],
    model: cfg?.model ?? '',
    temperature: options.temperature ?? 0.7,
    // DeepSeek 扩展参数：默认关闭思考模式（openai SDK 未内置该字段，需类型断言）
    thinking: { type: thinking ? 'enabled' : 'disabled' },
  } as OpenAI.ChatCompletionCreateParamsNonStreaming)

  const usage = completion.usage
  return {
    text: (completion.choices[0]?.message?.content ?? '').trim(),
    usage: usage
      ? {
          promptTokens: usage.prompt_tokens,
          completionTokens: usage.completion_tokens,
          totalTokens: usage.total_tokens,
        }
      : undefined,
  }
}
